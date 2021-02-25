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
use strum_macros::EnumIter;

/* crate use */
use crate::correct::exist::{Exist, Scenario};
use crate::set;

////////////////////////////////////
// Scenario for correct one error //
////////////////////////////////////
#[derive(Debug, EnumIter, Clone, Copy)]
pub enum ScenarioOne {
    I(usize, u8),
    S(usize, u8),
    D(usize, u8),
}

impl Scenario for ScenarioOne {
    fn init(&self, c: usize, k: u8) -> Self {
        match self {
            ScenarioOne::I(_, _) => ScenarioOne::I(c, k),
            ScenarioOne::S(_, _) => ScenarioOne::S(c, k),
            ScenarioOne::D(_, _) => ScenarioOne::D(c, k),
        }
    }

    fn c(&self) -> usize {
        match self {
            ScenarioOne::I(c, _) => *c,
            ScenarioOne::S(c, _) => *c,
            ScenarioOne::D(c, _) => *c,
        }
    }

    fn apply(&self, _valid_kmer: &set::BoxKmerSet, kmer: u64, _seq: &[u8]) -> Option<(u64, usize)> {
        match self {
            ScenarioOne::I(_, _) => Some((kmer, 2)),
            ScenarioOne::S(_, _) => Some((kmer, 1)),
            ScenarioOne::D(_, _) => Some((kmer, 0)),
        }
    }

    fn correct(&self, _valid_kmer: &set::BoxKmerSet, kmer: u64, _seq: &[u8]) -> (Vec<u8>, usize) {
        match self {
            ScenarioOne::I(_, _) => (vec![cocktail::kmer::bit2nuc(kmer & 0b11)], 2),
            ScenarioOne::S(_, _) => (vec![cocktail::kmer::bit2nuc(kmer & 0b11)], 1),
            ScenarioOne::D(_, _) => (vec![cocktail::kmer::bit2nuc(kmer & 0b11)], 0),
        }
    }
}

pub type One<'a> = Exist<'a, ScenarioOne>;

#[cfg(test)]
mod tests {

    use super::*;
    use crate::correct::Corrector;

    fn init() {
        let _ = env_logger::builder()
            .is_test(true)
            .filter_level(log::LevelFilter::Trace)
            .try_init();
    }

    fn filter<'a>(ori: &[u8]) -> Vec<u8> {
        ori.iter()
            .cloned()
            .filter(|x| *x != b'-')
            .collect::<Vec<u8>>()
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

        let refe = filter(b"ACTGA-CGAC");
        let read = filter(b"ACTGATCGAC");

        let mut data = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(&refe, 5) {
            data.set(kmer, true);
        }

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = One::new(&set, 2);

        assert_eq!(refe, corrector.correct(&read).as_slice());
        assert_eq!(refe, corrector.correct(&refe).as_slice());
    }

    #[test]
    fn cic_relaxe() {
        init();

        //                    GCGTAC-G
        //                   AGCGTAC
        //                      GTACTTG
        let refe = filter(b"GAGCGTAC-GTTGGAT");
        let read = filter(b"GAGCGTACTGTTGGAT");
        let conf = b"GCGTACGTGA";

        let mut data = pcon::solid::Solid::new(7);

        for kmer in cocktail::tokenizer::Tokenizer::new(&refe, 7) {
            data.set(kmer, true);
        }

        for kmer in cocktail::tokenizer::Tokenizer::new(conf, 7) {
            data.set(kmer, true);
        }

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = One::new(&set, 2);

        assert_eq!(refe, corrector.correct(&read).as_slice());
        assert_eq!(refe, corrector.correct(&refe).as_slice());
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
