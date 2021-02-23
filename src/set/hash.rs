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
pub use crate::set::KmerSet;

pub struct Hash {
    set: rustc_hash::FxHashSet<u64>,
    k: u8,
}

impl Hash {
    pub fn new<R>(inputs: Vec<R>, k: u8) -> Self
    where
        R: std::io::Read,
    {
        let mut set = rustc_hash::FxHashSet::default();

        for input in inputs {
            let mut records = bio::io::fasta::Reader::new(input).records();

            while let Some(Ok(record)) = records.next() {
                for cano in cocktail::tokenizer::Canonical::new(record.seq(), k) {
                    set.insert(cano);
                }
            }
        }

        Self { set, k }
    }
}

impl std::fmt::Debug for Hash {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Hash").finish()
    }
}

impl KmerSet for Hash {
    fn get(&self, kmer: u64) -> bool {
        self.set.contains(&cocktail::kmer::canonical(kmer, self.k))
    }

    fn k(&self) -> u8 {
        self.k
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    static FILE: &[u8] = b">1\nACGTGGGAATTGTGGCCACATCACGAGGTCCTGCGTATTGACGACTGTAAAGCGAGTGGCCGTGGAATTTCAAGCTCAATTAGCCGAACCAATCCGCCTA";

    #[test]
    fn canonical() {
        let file = std::io::Cursor::new(FILE);

        let hash = Hash::new(vec![file], 11);

        let set: crate::set::BoxKmerSet = Box::new(hash);

        let mut records = bio::io::fasta::Reader::new(FILE).records();
        for cano in cocktail::tokenizer::Canonical::new(records.next().unwrap().unwrap().seq(), 11)
        {
            assert!(set.get(cano))
        }
    }

    #[test]
    fn forward() {
        let file = std::io::Cursor::new(FILE);

        let hash = Hash::new(vec![file], 11);

        let set: crate::set::BoxKmerSet = Box::new(hash);

        let mut records = bio::io::fasta::Reader::new(FILE).records();
        for kmer in cocktail::tokenizer::Tokenizer::new(records.next().unwrap().unwrap().seq(), 11)
        {
            assert!(set.get(kmer))
        }
    }

    #[test]
    fn absence() {
        let file = std::io::Cursor::new(FILE);

        let hash = Hash::new(vec![file], 11);

        let set: crate::set::BoxKmerSet = Box::new(hash);

        assert!(!set.get(0));
    }
}
