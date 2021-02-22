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
use strum::IntoEnumIterator;
use strum_macros::EnumIter;

/* crate use */
use crate::correct::*;
use crate::set;

#[derive(Debug, EnumIter, Clone, Copy)]
enum Scenario {
    I,
    S,
    D,
}

impl Scenario {
    fn check(self) -> usize {
        match self {
            Scenario::I => 2,
            Scenario::S => 1,
            Scenario::D => 0,
        }
    }

    fn apply(self, kmer: u64, _seq: &[u8]) -> Option<(u64, usize)> {
        Some((kmer, self.check()))
    }

    fn correct(self, kmer: u64, _seq: &[u8]) -> (Vec<u8>, usize) {
        match self {
            Scenario::I => (vec![], 1),
            Scenario::S => (vec![cocktail::kmer::bit2nuc(kmer & 0b11)], 1),
            Scenario::D => (vec![cocktail::kmer::bit2nuc(kmer & 0b11)], 0),
        }
    }

    fn get_score(self, valid_kmer: &set::BoxKmerSet, ori: u64, seq: &[u8], c: usize) -> usize {
        if let Some((mut kmer, offset)) = self.apply(ori, seq) {
            if !valid_kmer.get(kmer) {
                return 0;
            }

            if offset + c > seq.len() {
                return 0;
            }

            let mut score = 0;

            for nuc in &seq[offset..offset + c] {
                kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(*nuc), valid_kmer.k());

                if valid_kmer.get(kmer) {
                    score += 1
                } else {
                    break;
                }
            }

            score
        } else {
            0
        }
    }
}

pub struct One<'a> {
    valid_kmer: &'a set::BoxKmerSet<'a>,
    c: u8,
}

impl<'a> One<'a> {
    pub fn new(valid_kmer: &'a set::BoxKmerSet, c: u8) -> Self {
        Self { valid_kmer, c }
    }

    fn get_scenarii(&self, kmer: u64, seq: &[u8]) -> Vec<Scenario> {
        let mut scenarii: Vec<Scenario> = Vec::new();

        for scenario in Scenario::iter() {
            if scenario.get_score(self.valid_kmer, kmer, seq, self.c as usize) == self.c as usize {
                scenarii.push(scenario)
            }
        }

        scenarii
    }
}

impl<'a> Corrector for One<'a> {
    fn valid_kmer(&self) -> &set::BoxKmerSet {
        self.valid_kmer
    }

    fn correct_error(&self, kmer: u64, seq: &[u8]) -> Option<(Vec<u8>, usize)> {
        let alts = alt_nucs(self.valid_kmer, kmer);

        if alts.len() != 1 {
            debug!("not one alts {:?}", alts);
            return None;
        }
        debug!("one alts {:?}", alts);

        let corr = add_nuc_to_end(kmer >> 2, alts[0], self.k());
        let mut scenarii = self.get_scenarii(corr, seq);

        if scenarii.is_empty() {
            debug!("no scenario");
            None
        } else if scenarii.len() == 1 {
            debug!("one scenario");
            Some(scenarii[0].correct(corr, seq))
        } else {
            scenarii.retain(|x| {
                last_is_valid(
                    self.valid_kmer,
                    corr,
                    &seq[(x.check() as usize)..(self.c as usize + x.check() + 1) as usize],
                )
            });

            if scenarii.len() == 1 {
                debug!("multi scenario one is better {:?}", scenarii[0]);
                Some(scenarii[0].correct(corr, seq))
            } else {
                debug!("multi scenario no better");
                None
            }
        }
    }
}

fn last_is_valid(valid_kmer: &set::BoxKmerSet, mut kmer: u64, nucs: &[u8]) -> bool {
    debug!("kmer {}", cocktail::kmer::kmer2seq(kmer, valid_kmer.k()));
    debug!("nucs {}", std::str::from_utf8(nucs).unwrap());
    nucs.iter()
        .for_each(|nuc| kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(*nuc), valid_kmer.k()));
    debug!(
        "last kmer {} valid {}",
        cocktail::kmer::kmer2seq(kmer, valid_kmer.k()),
        valid_kmer.get(kmer)
    );

    valid_kmer.get(kmer)
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

        let refe = b"GAGCGTACGTTGGAT";
        let read = b"GAGCGTACTGTTGGAT";
        let conf = b"GCGTACGTGA";

        let mut data = pcon::solid::Solid::new(7);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 7) {
            data.set(kmer, true);
        }

        for kmer in cocktail::tokenizer::Tokenizer::new(conf, 7) {
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

        let refe = b"ACTGACGACCC";
        let read = b"ACTGAGACCC";

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

        let refe = b"GAGCGTACGTTGGAT";
        let read = b"GAGCGTAGTTGGAT";
        let conf = b"GCGTACTT";

        let mut data = pcon::solid::Solid::new(7);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 7) {
            data.set(kmer, true);
        }

        for kmer in cocktail::tokenizer::Tokenizer::new(conf, 7) {
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
