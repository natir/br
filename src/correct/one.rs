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
use log::debug;

/* local use */
use crate::correct::*;
use crate::set;

pub struct One<'a> {
    valid_kmer: &'a set::BoxKmerSet<'a>,
    c: u8,
}

impl<'a> One<'a> {
    pub fn new(valid_kmer: &'a set::BoxKmerSet<'a>, c: u8) -> Self {
        Self { valid_kmer, c }
    }

    fn get_scenario(&self, kmer: u64, seq: &[u8]) -> Vec<(Vec<u8>, usize, usize)> {
        let mut scenario: Vec<(Vec<u8>, usize, usize)> = Vec::new();

        let alts = alt_nucs(self.valid_kmer, kmer);

        if alts.len() != 1 {
            debug!("not one alts {:?}", alts);
            return scenario;
        }
        debug!("one alts {:?}", alts);

        let corr = add_nuc_to_end(kmer >> 2, alts[0], self.k());

        if let Some(limit) = get_end_of_subseq(self.c as usize + 2, seq.len()) {
            // Substitution
	    let score = get_kmer_score(self.valid_kmer, corr, &seq[1..limit]);
            if score >= self.c as usize {
                debug!("it's a substitution {}", alts[0]);
                scenario.push((vec![cocktail::kmer::bit2nuc(alts[0])], 1, score));
            }
        }

        if let Some(limit) = get_end_of_subseq(self.c as usize + 3, seq.len()) {
            // Insertion
	    let score = get_kmer_score(self.valid_kmer, corr, &seq[2..limit]);
            if  score >= self.c as usize {
                debug!("it's a insertion");
                scenario.push((Vec::new(), 1, score));
            }
        }

        if let Some(limit) = get_end_of_subseq(self.c as usize + 1, seq.len()) {
            // Deletion
	    let score = get_kmer_score(self.valid_kmer, corr, &seq[0..limit]);
            if  score >= self.c as usize {
                debug!("it's a deletion {}", alts[0]);
                scenario.push((vec![cocktail::kmer::bit2nuc(alts[0])], 0, score));
            }
        }

        scenario
    }
}

impl<'a> Corrector for One<'a> {
    fn valid_kmer(&self) -> &set::BoxKmerSet {
        self.valid_kmer
    }

    fn correct_error(&self, kmer: u64, seq: &[u8]) -> Option<(Vec<u8>, usize)> {
        let mut scenario = self.get_scenario(kmer, seq);
	
        if scenario.is_empty() {
            debug!("no scenario {:?}", scenario);
            None
        } else if scenario.len() == 1 {
            debug!("one scenario");
            Some((scenario[0].0.clone(), scenario[0].1))
        } else {
	    scenario.retain(|x| x.2 == (self.c as usize + 1));
	    if scenario.len() == 1 {
		debug!("multi scenario one good");
		Some((scenario[0].0.clone(), scenario[0].1))
	    } else {
		debug!("multi scenario no good");
		None
	    }
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

fn get_kmer_score(valid_kmer: &set::BoxKmerSet, mut kmer: u64, nucs: &[u8]) -> usize {
    let mut score = 0;

    for nuc in nucs {
        kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(*nuc), valid_kmer.k());

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

        let refe = b"ACTGACGAC";
        let read = b"ACTGATGAC";

        let mut data = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = One::new(&set, 2);

        assert_eq!(refe, corrector.correct(read).as_slice());
        assert_eq!(refe, corrector.correct(refe).as_slice());
    }

    #[test]
    fn csc_relaxe() {
        init();

        let refe = b"ACTGACCACT";
        let read = b"ACTGATCACT";
        let conf = b"ACTGACAC";

        let mut data = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        for kmer in cocktail::tokenizer::Tokenizer::new(conf, 5) {
            data.set(kmer, true);
        }

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = One::new(&set, 2);

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

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = One::new(&set, 2);

        assert_eq!(read, corrector.correct(read).as_slice()); // don't correct
        assert_eq!(refe, corrector.correct(refe).as_slice());
    }

    #[test]
    fn cic() {
        init();

        let refe = b"ACTGACGAC";
        let read = b"ACTGATCGAC";

        let mut data = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = One::new(&set, 2);

        assert_eq!(refe, corrector.correct(read).as_slice());
        assert_eq!(refe, corrector.correct(refe).as_slice());
    }

    #[test]
    fn cic_relaxe() {
        init();

        let refe = b"ACTGACGAC";
        let read = b"ACTGATCGAC";
        let conf = b"ACTGACGAC";

        let mut data = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        for kmer in cocktail::tokenizer::Tokenizer::new(conf, 5) {
            data.set(kmer, true);
        }

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = One::new(&set, 2);

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

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = One::new(&set, 2);

        assert_eq!(read, corrector.correct(read).as_slice()); // don't correct
        assert_eq!(refe, corrector.correct(refe).as_slice());
    }

    #[test]
    fn cdc() {
        init();

        let refe = b"ACTGACGAC";
        let read = b"ACTGAGAC";

        let mut data = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = One::new(&set, 2);

        assert_eq!(refe, corrector.correct(read).as_slice());
        assert_eq!(refe, corrector.correct(refe).as_slice());
    }

    #[test]
    fn cdc_relaxe() {
        init();

        let refe = b"ACTGACGAC";
        let read = b"ACTGAGAC";
        let conf = b"ACTGACTAC";

        let mut data = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        for kmer in cocktail::tokenizer::Tokenizer::new(conf, 5) {
            data.set(kmer, true);
        }

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = One::new(&set, 2);

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

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = One::new(&set, 2);

        assert_eq!(read, corrector.correct(read).as_slice()); // don't correct
        assert_eq!(refe, corrector.correct(refe).as_slice());
    }
}
