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

/* local use */
use crate::set;

const MASK_LOOKUP: [u64; 32] = {
    let mut lookup = [0; 32];

    let mut k = 1;
    while k < 32 {
        lookup[k] = (1 << (2 * k)) - 1;

        k += 1;
    }

    lookup
};

#[inline(always)]
pub(crate) fn mask(k: u8) -> u64 {
    MASK_LOOKUP[k as usize]
}

pub trait Corrector {
    fn valid_kmer(&self) -> &set::BoxKmerSet;

    fn correct_error(&self, kmer: u64, seq: &[u8]) -> Option<(Vec<u8>, usize)>;

    fn k(&self) -> u8 {
        self.valid_kmer().k()
    }

    fn correct(&self, seq: &[u8]) -> Vec<u8> {
        let mut correct: Vec<u8> = Vec::with_capacity(seq.len());

        if seq.len() < self.k() as usize {
            return seq.to_vec();
        }

        let mut i = self.k() as usize;
        let mut kmer = cocktail::kmer::seq2bit(&seq[0..i]);

        for n in &seq[0..i] {
            correct.push(*n);
        }

        let mut previous = self.valid_kmer().get(kmer);
        while i < seq.len() {
            let nuc = seq[i];

            kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(nuc), self.k());

            if !self.valid_kmer().get(kmer) && previous {
                if let Some((local_correct, offset)) = self.correct_error(kmer, &seq[i..]) {
                    kmer >>= 2;

                    for nuc in local_correct {
                        kmer = add_nuc_to_end(
                            kmer,
                            cocktail::kmer::nuc2bit(nuc),
                            self.valid_kmer().k(),
                        );
                        correct.push(nuc);
                    }

                    log::debug!("error at position {} cor", i);

                    previous = true;
                    i += offset;
                } else {
                    correct.push(nuc);

                    log::debug!("error at position {} not", i);

                    i += 1;
                    previous = false;
                }
            } else {
                previous = self.valid_kmer().get(kmer);
                correct.push(nuc);

                i += 1;
            }
        }

        correct
    }
}

pub(crate) fn add_nuc_to_end(kmer: u64, nuc: u64, k: u8) -> u64 {
    ((kmer << 2) & mask(k)) ^ nuc
}

pub(crate) fn alt_nucs(valid_kmer: &set::BoxKmerSet, ori: u64) -> Vec<u64> {
    next_nucs(valid_kmer, ori >> 2)
}

pub(crate) fn next_nucs(valid_kmer: &set::BoxKmerSet, kmer: u64) -> Vec<u64> {
    let mut correct_nuc: Vec<u64> = Vec::with_capacity(4);

    for alt_nuc in 0..4 {
        if valid_kmer.get(add_nuc_to_end(kmer, alt_nuc, valid_kmer.k())) {
            correct_nuc.push(alt_nuc);
        }
    }

    correct_nuc
}

pub(crate) fn error_len(
    subseq: &[u8],
    mut kmer: u64,
    valid_kmer: &set::BoxKmerSet,
) -> (usize, u64) {
    let mut j = 0;

    loop {
        j += 1;

        if j >= subseq.len() {
            break;
        }

        kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(subseq[j]), valid_kmer.k());

        if valid_kmer.get(kmer) {
            break;
        }
    }

    (j, kmer)
}

pub mod exist;
pub mod gap_size;
pub mod graph;
pub mod greedy;

pub use exist::one::One;
pub use exist::two::Two;
pub use gap_size::GapSize;
pub use graph::Graph;
pub use greedy::Greedy;

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn found_alt_kmer() {
        let mut data = pcon::solid::Solid::new(5);
        data.set(cocktail::kmer::seq2bit(b"ACTGA"), true);
        data.set(cocktail::kmer::seq2bit(b"ACTGT"), true);

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let kmer = cocktail::kmer::seq2bit(b"ACTGC");

        assert_eq!(alt_nucs(&set, kmer), vec![0, 2]);
    }
}
