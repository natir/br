/*
Copyright (c) 2020 Pierre Marijon <pierre.marijon@hhu.de>

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
use crate::set::KmerSet;

pub struct Pcon {
    set: pcon::solid::Solid,
}

impl Pcon {
    pub fn new(set: pcon::solid::Solid) -> Self {
        Pcon { set }
    }
}

impl KmerSet for Pcon {
    fn get(&self, kmer: u64) -> bool {
        self.set.get(kmer)
    }

    fn k(&self) -> u8 {
        self.set.k
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    static SEQ: &[u8] = b"ACGTGGGAATTGTGGCCACATCACGAGGTCCTGCGTATTGACGACTGTAAAGCGAGTGGCCGTGGAATTTCAAGCTCAATTAGCCGAACCAATCCGCCTA";

    #[test]
    fn canonical() {
        let mut solid = pcon::solid::Solid::new(11);
        for cano in cocktail::tokenizer::Canonical::new(SEQ, 11) {
            solid.set(cano, true);
        }

        let set: crate::set::BoxKmerSet = Box::new(Pcon::new(solid));

        for cano in cocktail::tokenizer::Canonical::new(SEQ, 11) {
            assert!(set.get(cano))
        }
    }

    #[test]
    fn forward() {
        let mut solid = pcon::solid::Solid::new(11);
        for cano in cocktail::tokenizer::Canonical::new(SEQ, 11) {
            solid.set(cano, true);
        }

        let set: crate::set::BoxKmerSet = Box::new(Pcon::new(solid));

        for kmer in cocktail::tokenizer::Tokenizer::new(SEQ, 11) {
            assert!(set.get(kmer))
        }
    }

    #[test]
    fn absence() {
        let mut solid = pcon::solid::Solid::new(11);
        for cano in cocktail::tokenizer::Canonical::new(SEQ, 11) {
            solid.set(cano, true);
        }

        let set: crate::set::BoxKmerSet = Box::new(Pcon::new(solid));

        assert!(!set.get(0));
    }
}
