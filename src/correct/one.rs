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
use log::{debug, info};

/* local use */
use crate::correct::*;

pub struct One<'a> {
    valid_kmer: &'a pcon::solid::Solid,
    c: u8,
}

impl<'a> One<'a> {
    pub fn new(valid_kmer: &'a pcon::solid::Solid, c: u8) -> Self {
        Self { valid_kmer, c }
    }
}

impl<'a> Corrector for One<'a> {
    fn valid_kmer(&self) -> &pcon::solid::Solid {
	self.valid_kmer
    }
    
    fn correct_error(&self, mut kmer: u64, seq: &[u8]) -> Option<(Vec<u8>, usize)> {
	let alts = alt_nucs(self.valid_kmer(), kmer);

	if alts.len() != 1 {
            debug!("not one alts {:?}", alts);
            return None;
	}
	debug!("one alts {:?}", alts);

	let corr = add_nuc_to_end(kmer >> 2, alts[0], self.k());

	let mut scenario: Vec<(Vec<u8>, usize)> = Vec::new();

	if let Some(limit) = get_end_of_subseq(self.c as usize + 1, seq.len()) {
            // Substitution
            if get_kmer_score(self.valid_kmer(), corr, &seq[1..limit]) == self.c {
		debug!("it's a substitution {}", alts[0]);
		scenario.push((vec![cocktail::kmer::bit2nuc(alts[0])], 1));
            }
	}

	if let Some(limit) = get_end_of_subseq(self.c as usize + 2, seq.len()) {
            // Insertion
            if get_kmer_score(self.valid_kmer(), corr, &seq[2..limit]) == self.c {
		debug!("it's a insertion");
		scenario.push((Vec::new(), 1));
            }
	}

	if let Some(limit) = get_end_of_subseq(self.c as usize, seq.len()) {
            // Deletion
            if get_kmer_score(self.valid_kmer(), corr, &seq[0..limit]) == self.c {
		debug!("it's a deletion {}", alts[0]);
		scenario.push((vec![cocktail::kmer::bit2nuc(alts[0])], 0));
            }
	}

	if scenario.len() == 1 {
            debug!("one scenario");
            scenario.pop()
	} else {
            debug!("multi scenario {:?}", scenario);
            None
	}

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
        kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(*nuc), valid_kmer.k);

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

        let corrector = One::new(&data, 2);

        assert_eq!(refe, corrector.correct(read).as_slice());
        assert_eq!(refe, corrector.correct(refe).as_slice());
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

        let corrector = One::new(&data, 2);

        assert_eq!(read, corrector.correct(read).as_slice()); // don't correct
        assert_eq!(refe, corrector.correct(refe).as_slice());
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

        let corrector = One::new(&data, 2);

        assert_eq!(refe, corrector.correct(read).as_slice());
        assert_eq!(refe, corrector.correct(refe).as_slice());
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

        let corrector = One::new(&data, 2);

        assert_eq!(read, corrector.correct(read).as_slice()); // don't correct
        assert_eq!(refe, corrector.correct(refe).as_slice());
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

        let corrector = One::new(&data, 2);

        assert_eq!(refe, corrector.correct(read).as_slice());
        assert_eq!(refe, corrector.correct(refe).as_slice());
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

        let corrector = One::new(&data, 2);

        assert_eq!(read, corrector.correct(read).as_slice()); // don't correct
        assert_eq!(refe, corrector.correct(refe).as_slice());
    }
}
