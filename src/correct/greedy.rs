/*
Copyright (c) 2020 Pierre Marijon <pmarijon@mmci.uni-saarland.de>

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

pub fn correct(seq: &[u8], valid_kmer: &pcon::solid::Solid, c: u8) -> Vec<u8> {
    let mut correct: Vec<u8> = Vec::with_capacity(seq.len());

    if seq.len() < valid_kmer.k as usize {
        return seq.to_vec();
    }

    let mut i = valid_kmer.k as usize;
    let mut kmer = cocktail::kmer::seq2bit(&seq[0..i]);

    for n in &seq[0..i] {
        correct.push(*n);
    }

    let mut previous = true;
    while i < seq.len() {
        let nuc = seq[i];
	trace!("nucleotide {}", nuc);

        kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(nuc));

        if !valid_kmer.get(kmer) && previous {
            debug!("kmer {} isn't exist", cocktail::kmer::kmer2seq(kmer, valid_kmer.k));

            if let Some((corr, offset)) = correct_error(kmer, &seq[i..], c, valid_kmer) {
                previous = true;

		kmer >>= 2;
                for nuc in corr {
                    kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(nuc));

                    correct.push(nuc);
                }

                i += offset;

                debug!("error at position {} has been corrected", i);
            } else {
                previous = false;
                correct.push(nuc);

                i += 1;

                debug!("error at position {} hasn't been corrected", i);
            }
        } else {
            previous = valid_kmer.get(kmer);
            correct.push(nuc);

            i += 1;
        }
    }

    correct
}

fn correct_error(
    kmer: u64,
    seq: &[u8],
    th: u8,
    valid_kmer: &pcon::solid::Solid,
) -> Option<(Vec<u8>, usize)> {
    let succs = alt_nucs(valid_kmer, kmer);

    if succs.len() != 1 {
        return None;
    }

    let corr = add_nuc_to_end(kmer >> 2, succs[0]);

    let mut scenario: Vec<(Vec<u8>, usize)> = Vec::new();

    if let Some(limit) = get_end_of_subseq(th as usize + 1, seq.len()) {
        // Substitution
        if get_kmer_score(valid_kmer, corr, &seq[1..limit]) == th {
	    trace!("it's a substitution {}", succs[0]);
            scenario.push((vec![cocktail::kmer::bit2nuc(succs[0])], 1));
        }
    }

    if let Some(limit) = get_end_of_subseq(th as usize + 2, seq.len()) {
        // Insertion
        if get_kmer_score(valid_kmer, corr, &seq[2..limit]) == th {
	    trace!("it's a insertion");
            scenario.push((Vec::new(), 1));
        }
    }

    if let Some(limit) = get_end_of_subseq(th as usize, seq.len()) {
        // Deletion
        if get_kmer_score(valid_kmer, corr, &seq[0..limit]) == th {
	    trace!("it's a deletion {}", succs[0]);
            scenario.push((vec![cocktail::kmer::bit2nuc(succs[0])], 0));
        }
    }

    if scenario.len() == 1 {
	trace!("one scenario");
        scenario.pop()
    } else {
	trace!("multi scenario {}", scenario.len());
        None
    }
}

fn get_end_of_subseq(offset: usize, max_length: usize) -> Option<usize> {
    if offset > max_length {
	None	    
    } else {
        Some(offset)
    }
}

fn get_kmer_score(valid_kmer: &pcon::solid::Solid, mut kmer: u64, nucs: &[u8]) -> u8 {
    let mut score = 0;

    for nuc in nucs {
        kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(*nuc));

        if valid_kmer.get(kmer) {
            score += 1
        }
    }

    score
}

#[cfg(test)]
mod tests {

    use super::*;

    fn init() {
        let _ = env_logger::builder().is_test(true).try_init();

        init_masks(5);
    }

    #[test]
    fn csc() {
        init();

        let refe = b"ACTGACGA";
        let read = b"ACTGATGA";

        let mut data = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        assert_eq!(refe, correct(read, &data, 2).as_slice());
        assert_eq!(refe, correct(refe, &data, 2).as_slice());
    }

    #[test]
    fn cssc() {
        init();

        let refe = b"ACTGACGAG";
        let read = b"ACTGATAAG";

        let mut data = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        assert_eq!(read, correct(read, &data, 2).as_slice()); // don't correct
        assert_eq!(refe, correct(refe, &data, 2).as_slice());
    }

    
    #[test]
    fn cic() {
        init();

        let refe = b"ACTGACGA";
        let read = b"ACTGATCGA";

        let mut data = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        assert_eq!(refe, correct(read, &data, 2).as_slice());
        assert_eq!(refe, correct(refe, &data, 2).as_slice());
    }

    #[test]
    fn ciic() {
        init();

        let refe = b"ACTGACGA";
        let read = b"ACTGATTCGA";

        let mut data = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        assert_eq!(read, correct(read, &data, 2).as_slice()); // don't correct
        assert_eq!(refe, correct(refe, &data, 2).as_slice());
    }
    
    #[test]
    fn cdc() {
        init();

        let refe = b"ACTGACGA";
        let read = b"ACTGAGA";

        let mut data = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        assert_eq!(refe, correct(read, &data, 2).as_slice());
        assert_eq!(refe, correct(refe, &data, 2).as_slice());
    }

    #[test]
    fn cddc() {
        init();

        let refe = b"ACTGACGAG";
        let read = b"ACTGAAG";

        let mut data = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        assert_eq!(read, correct(read, &data, 2).as_slice()); // don't correct
        assert_eq!(refe, correct(refe, &data, 2).as_slice());
    }
}
