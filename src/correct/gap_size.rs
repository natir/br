/*
Copyright (c) 2020 Pierre Marijon <pmarijon@hhu.de>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */

/* crate use */
use log::{debug, info};

/* local use */
use crate::correct::*;

pub fn correct(seq: &[u8], valid_kmer: &pcon::solid::Solid, c: u8) -> Vec<u8> {
    let mut correct;
    let mut i;
    let mut kmer;

    if let Some(res) = init_correction(seq, valid_kmer.k) {
        correct = res.0;
        i = res.1;
        kmer = res.2;
    } else {
        return seq.to_vec();
    }

    let mut previous = valid_kmer.get(kmer);
    while i < seq.len() {
        let nuc = seq[i];

        kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(nuc));

        if previous && !valid_kmer.get(kmer) {
            let (error_len, first_correct_kmer) = error_len(&seq[i..], kmer, valid_kmer);

            if error_len < valid_kmer.k as usize {
                // Delettion
                if let Some(local_correct) =
                    graph::correct_error(kmer, first_correct_kmer, valid_kmer)
                {
                    correct.extend(local_correct.iter());
                    previous = true;

                    info!(
                        "error at position {} of length {} deletion cor",
                        i, error_len
                    );
                } else {
                    correct.extend(&seq[i..i + error_len]);
                    previous = false;
                    info!(
                        "error at position {} of length {} deletion not",
                        i, error_len
                    );
                }
                kmer = first_correct_kmer;
                i += error_len + 1;
            } else if error_len > valid_kmer.k as usize {
                // insertion substitution large

                if let Some(local_corr) =
                    correct_error(kmer, valid_kmer, error_len - valid_kmer.k as usize)
                {
                    info!(
                        "error at position {} of length {} sub/ins large cor",
                        i, error_len
                    );

                    kmer >>= 2;

                    for nuc in local_corr.iter() {
                        correct.push(*nuc);
                        kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(*nuc));
                        i += 1;
                    }

                    // TODO: side effect can't detect insertion at end of reads check it
                    if i + local_corr.len() <= seq[i..].len() {
                        if local_corr == &seq[i..(i + local_corr.len())] {
                            debug!("It's an insertion");

                            i += local_corr.len();
                        }
                    }
                    previous = true;
                } else {
                    for nuc in &seq[i..(i + error_len)] {
                        correct.push(*nuc);
                        kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(*nuc));
                        i += 1;
                    }
                    previous = false;

                    info!(
                        "error at position {} of length {} sub/ins large not",
                        i, error_len
                    );
                }
            } else if error_len == valid_kmer.k as usize {
                // insertion substitution 1

                if let Some((corr, offset)) = greedy::correct_error(kmer, &seq[i..], c, valid_kmer)
                {
                    kmer >>= 2;
                    for nuc in corr {
                        kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(nuc));

                        correct.push(nuc);
                    }

                    i += offset;

                    previous = true;

                    info!(
                        "error at position {} of length {} ins/sub 1 cor",
                        i, error_len
                    );
                } else {
                    correct.push(nuc);

                    i += 1;

                    previous = false;

                    info!(
                        "error at position {} of length {} ins/sub 1 not",
                        i, error_len
                    );
                }
            }
        } else {
            previous = valid_kmer.get(kmer);

            correct.push(nuc);

            i += 1;
        }
    }

    correct
}

pub(crate) fn correct_error(
    kmer: u64,
    valid_kmer: &pcon::solid::Solid,
    gap_size: usize,
) -> Option<Vec<u8>> {
    let mut alts = alt_nucs(valid_kmer, kmer);

    if alts.len() != 1 {
        debug!("not one alts {:?}", alts);
        return None;
    }

    let mut corr = add_nuc_to_end(kmer >> 2, alts[0]);
    let mut local_corr = vec![cocktail::kmer::bit2nuc(alts[0])];

    for _ in 0..gap_size {
        debug!("kmer {:?}", cocktail::kmer::kmer2seq(corr, valid_kmer.k));

        alts = next_nucs(valid_kmer, corr);

        if alts.len() != 1 {
            debug!("failled multiple successor {:?}", alts);
            return None;
        }

        corr = add_nuc_to_end(corr, alts[0]);

        local_corr.push(cocktail::kmer::bit2nuc(alts[0]));
    }

    Some(local_corr)
}

#[cfg(test)]
mod tests {

    use super::*;

    fn init() {
        let _ = env_logger::builder()
            .is_test(true)
            .filter_level(log::LevelFilter::Trace)
            .try_init();

        init_masks(5);
    }

    #[test]
    fn found_alt_kmer() {
        init();

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        data.set(cocktail::kmer::seq2bit(b"ACTGA"), true);
        data.set(cocktail::kmer::seq2bit(b"ACTGT"), true);

        let kmer = cocktail::kmer::seq2bit(b"ACTGC");

        assert_eq!(alt_nucs(&data, kmer), vec![0, 2]);
    }

    #[test]
    fn csc() {
        init();

        let refe = b"AGCGTATCTT";
        //           ||||| |||||
        let read = b"AGCGTTTCTT";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        assert_eq!(refe, correct(read, &data, 2).as_slice()); // test correction work
        assert_eq!(refe, correct(refe, &data, 2).as_slice()); // test not overcorrection
    }

    #[test]
    fn cssc() {
        init();

        let refe = b"TCTCTAATCTTC";
        //           |||||  |||||
        let read = b"TCTCTGGTCTTC";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        assert_eq!(refe, correct(read, &data, 2).as_slice()); // test correction work
        assert_eq!(refe, correct(refe, &data, 2).as_slice()); // test not overcorrection
    }

    #[test]
    fn csssc() {
        init();

        let refe = b"TCTCTAAATCTTC";
        //           |||||  |||||
        let read = b"TCTCTGGGTCTTC";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        assert_eq!(refe, correct(read, &data, 2).as_slice()); // test correction work
        assert_eq!(refe, correct(refe, &data, 2).as_slice()); // test not overcorrect
    }

    #[test]
    fn cscsc() {
        init();

        let refe = b"TCTTTACATTTTT";
        //           ||||| | |||||
        let read = b"TCTTTGCGTTTTT";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        assert_eq!(refe, correct(read, &data, 2).as_slice()); // test correction work
        assert_eq!(refe, correct(refe, &data, 2).as_slice()); // test not overcorrection
    }

    #[test]
    fn cdc() {
        init();

        let refe = b"GATACATGGACACTAGTATG";
        //           ||||||||||
        let read = b"GATACATGGAACTAGTATG";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        assert_eq!(refe, correct(read, &data, 2).as_slice());
        assert_eq!(refe, correct(refe, &data, 2).as_slice());
    }

    #[test]
    fn cddc() {
        init();

        let refe = b"CAAAGCATTTTT";
        //           |||||
        let read = b"CAAAGTTTTT";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        assert_eq!(refe, correct(read, &data, 2).as_slice());
        assert_eq!(refe, correct(refe, &data, 2).as_slice());
    }

    #[test]
    fn cic() {
        init();

        let refe = b"GGATAACTCT";
        //           |||||
        let read = b"GGATATACTCT";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        assert_eq!(refe, correct(read, &data, 2).as_slice()); // test correction work
        assert_eq!(refe, correct(refe, &data, 2).as_slice()); // test not overcorrection
    }

    #[test]
    fn ciic() {
        init();

        let refe = b"GATACATGGACACTAGTATGACGA";
        //           ||||||||||
        let read = b"GATACATGGATTCACTAGTATGACGA";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        assert_eq!(refe, correct(read, &data, 2).as_slice()); // test correction work
        assert_eq!(refe, correct(refe, &data, 2).as_slice()); // test not overcorrection
    }
}