/*
Copyright (c) 2019 Pierre Marijon <pmarijon@mmci.uni-saarland.de>

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
use log::{debug, trace};

/* local use */
use crate::correct::*;

pub fn correct(seq: &[u8], valid_kmer: &pcon::solid::Solid) -> Vec<u8> {
    let mut correct: Vec<u8> = Vec::with_capacity(seq.len());

    if seq.len() < valid_kmer.k as usize {
        return seq.to_vec();
    }

    let mut i = valid_kmer.k as usize;
    let mut kmer = cocktail::kmer::seq2bit(&seq[0..i]);

    for n in &seq[0..i] {
        correct.push(*n);
    }

    while i < seq.len() {
        let nuc = seq[i];

        kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(nuc));

        if !valid_kmer.get(kmer) {
            trace!("kmer {} isn't exist", kmer);

            let (error_len, first_correct_kmer) = error_len(&seq[i..], kmer, valid_kmer);

            if let Some(local_correct) = correct_error(kmer, first_correct_kmer, valid_kmer) {
                correct.extend(local_correct.iter());
                trace!(
                    "error at position {} of length {} has been corrected",
                    i,
                    error_len
                );
            } else {
                correct.extend(&seq[i..i + error_len]);
                trace!(
                    "error at position {} of length {} hasn't been corrected",
                    i,
                    error_len
                );
            }

            kmer = first_correct_kmer;
            i += error_len + 1;
        } else {
            correct.push(nuc);

            i += 1;
        }
    }

    correct
}

fn correct_error(
    mut kmer: u64,
    first_correct_kmer: u64,
    valid_kmer: &pcon::solid::Solid,
) -> Option<Vec<u8>> {
    let mut local_corr = Vec::new();

    let succs = alt_nucs(valid_kmer, kmer);
    if succs.len() != 1 {
        return None;
    }

    kmer = add_nuc_to_end(kmer >> 2, succs[0]);
    local_corr.push(cocktail::kmer::bit2nuc(succs[0]));

    while valid_kmer.get(kmer) {
        let succs = next_nucs(valid_kmer, kmer);

        if succs.len() != 1 {
            return None;
        }

        kmer = add_nuc_to_end(kmer, succs[0]);

        local_corr.push(cocktail::kmer::bit2nuc(succs[0]));

        if kmer == first_correct_kmer {
            break;
        }
    }

    Some(local_corr)
}

fn error_len(subseq: &[u8], mut kmer: u64, valid_kmer: &pcon::solid::Solid) -> (usize, u64) {
    let mut j = 0;

    loop {
        j += 1;

        if j >= subseq.len() {
            break;
        }

        kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(subseq[j]));

        if valid_kmer.get(kmer) {
            break;
        }
    }

    (j, kmer)
}

#[cfg(test)]
mod tests {

    use super::*;

    fn init() {
        let _ = env_logger::builder().is_test(true).try_init();

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

        let refe = b"TCTTTATTTTC";
        //           ||||| |||||
        let read = b"TCTTTGTTTTC";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        assert_eq!(refe, correct(read, &data).as_slice()); // test correction work
        assert_eq!(refe, correct(refe, &data).as_slice()); // test not overcorrection
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

        assert_eq!(refe, correct(read, &data).as_slice()); // test correction work
        assert_eq!(refe, correct(refe, &data).as_slice()); // test not overcorrection
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

        assert_eq!(refe, correct(read, &data).as_slice()); // test correction work
        assert_eq!(refe, correct(refe, &data).as_slice()); // test not overcorrect
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

        assert_eq!(refe, correct(read, &data).as_slice()); // test correction work
        assert_eq!(refe, correct(refe, &data).as_slice()); // test not overcorrection
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

        assert_eq!(refe, correct(read, &data).as_slice());
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

        assert_eq!(refe, correct(read, &data).as_slice());
    }

    #[test]
    fn cic() {
        init();

        let refe = b"GATACATGGACACTAGTATG";
        //           ||||||||||
        let read = b"GATACATGGATCACTAGTATG";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        assert_eq!(refe, correct(read, &data).as_slice()); // test correction work
        assert_eq!(refe, correct(refe, &data).as_slice()); // test not overcorrection
    }

    #[test]
    fn ciic() {
        init();

        let refe = b"GATACATGGACACTAGTATG";
        //           ||||||||||
        let read = b"GATACATGGATTCACTAGTATG";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        assert_eq!(refe, correct(read, &data).as_slice()); // test correction work
        assert_eq!(refe, correct(refe, &data).as_slice()); // test not overcorrection
    }
}
