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

/* crate use */
use strum_macros::EnumIter;

/* crate use */
use crate::correct::exist::{Exist, Scenario};
use crate::correct::*;
use crate::set;

////////////////////////////////////
// Scenario for correct two error //
////////////////////////////////////
#[derive(Debug, EnumIter, Clone, Copy)]
pub enum ScenarioTwo {
    I(usize, u8),
    S(usize, u8),
    D(usize, u8),

    II(usize, u8),
    IS(usize, u8),
    SS(usize, u8),
    SD(usize, u8),
    DD(usize, u8),

    ICI(usize, u8),
    ICS(usize, u8),
    ICD(usize, u8),
    SCI(usize, u8),
    SCS(usize, u8),
    SCD(usize, u8),
    DCI(usize, u8),
    DCS(usize, u8),
    DCD(usize, u8),
}

impl Scenario for ScenarioTwo {
    fn init(&self, c: usize, k: u8) -> Self {
        match self {
            ScenarioTwo::I(_, _) => ScenarioTwo::I(c, k),
            ScenarioTwo::S(_, _) => ScenarioTwo::S(c, k),
            ScenarioTwo::D(_, _) => ScenarioTwo::D(c, k),
            ScenarioTwo::II(_, _) => ScenarioTwo::II(c, k),
            ScenarioTwo::IS(_, _) => ScenarioTwo::IS(c, k),
            ScenarioTwo::SS(_, _) => ScenarioTwo::SS(c, k),
            ScenarioTwo::SD(_, _) => ScenarioTwo::SD(c, k),
            ScenarioTwo::DD(_, _) => ScenarioTwo::DD(c, k),
            ScenarioTwo::ICI(_, _) => ScenarioTwo::ICI(c, k),
            ScenarioTwo::ICS(_, _) => ScenarioTwo::ICS(c, k),
            ScenarioTwo::ICD(_, _) => ScenarioTwo::ICD(c, k),
            ScenarioTwo::SCI(_, _) => ScenarioTwo::SCI(c, k),
            ScenarioTwo::SCS(_, _) => ScenarioTwo::SCS(c, k),
            ScenarioTwo::SCD(_, _) => ScenarioTwo::SCD(c, k),
            ScenarioTwo::DCI(_, _) => ScenarioTwo::DCI(c, k),
            ScenarioTwo::DCS(_, _) => ScenarioTwo::DCS(c, k),
            ScenarioTwo::DCD(_, _) => ScenarioTwo::DCD(c, k),
        }
    }

    fn c(&self) -> usize {
        match self {
            ScenarioTwo::I(c, _) => *c,
            ScenarioTwo::S(c, _) => *c,
            ScenarioTwo::D(c, _) => *c,
            ScenarioTwo::II(c, _) => *c,
            ScenarioTwo::IS(c, _) => *c,
            ScenarioTwo::SS(c, _) => *c,
            ScenarioTwo::SD(c, _) => *c,
            ScenarioTwo::DD(c, _) => *c,
            ScenarioTwo::ICI(c, _) => *c,
            ScenarioTwo::ICS(c, _) => *c,
            ScenarioTwo::ICD(c, _) => *c,
            ScenarioTwo::SCI(c, _) => *c,
            ScenarioTwo::SCS(c, _) => *c,
            ScenarioTwo::SCD(c, _) => *c,
            ScenarioTwo::DCI(c, _) => *c,
            ScenarioTwo::DCS(c, _) => *c,
            ScenarioTwo::DCD(c, _) => *c,
        }
    }

    fn apply(
        &self,
        valid_kmer: &set::BoxKmerSet,
        mut kmer: u64,
        seq: &[u8],
    ) -> Option<(u64, usize)> {
        match self {
            ScenarioTwo::I(_, _) => Some((kmer, 2)),
            ScenarioTwo::S(_, _) => Some((kmer, 1)),
            ScenarioTwo::D(_, _) => Some((kmer, 0)),
            ScenarioTwo::II(_, _) => Some((kmer, 3)),
            ScenarioTwo::IS(_, k) => {
                if seq.len() < 2 {
                    None
                } else {
                    kmer = add_nuc_to_end(kmer >> 2, cocktail::kmer::nuc2bit(seq[1]), *k);

                    if valid_kmer.get(kmer) {
                        None
                    } else {
                        let alts = alt_nucs(valid_kmer, kmer);

                        if alts.len() != 1 {
                            None
                        } else {
                            Some((add_nuc_to_end(kmer >> 2, alts[0], *k), 3))
                        }
                    }
                }
            }
            ScenarioTwo::SS(_, _) => None,
            ScenarioTwo::SD(_, _) => None,
            ScenarioTwo::DD(_, _) => None,
            ScenarioTwo::ICI(_, _) => None,
            ScenarioTwo::ICS(_, _) => None,
            ScenarioTwo::ICD(_, _) => None,
            ScenarioTwo::SCI(_, _) => None,
            ScenarioTwo::SCS(_, _) => None,
            ScenarioTwo::SCD(_, _) => None,
            ScenarioTwo::DCI(_, _) => None,
            ScenarioTwo::DCS(_, _) => None,
            ScenarioTwo::DCD(_, _) => None,
        }
    }

    fn correct(&self, kmer: u64, _seq: &[u8]) -> (Vec<u8>, usize) {
        match self {
            ScenarioTwo::I(_, _) => (vec![], 2),
            ScenarioTwo::S(_, _) => (vec![cocktail::kmer::bit2nuc(kmer & 0b11)], 1),
            ScenarioTwo::D(_, _) => (vec![cocktail::kmer::bit2nuc(kmer & 0b11)], 0),
            ScenarioTwo::II(_, _) => (vec![], 2),
            ScenarioTwo::IS(_, _) => (vec![cocktail::kmer::bit2nuc(kmer & 0b11)], 2),

            /*    ScenarioTwo::SS => {
                None
                },
                ScenarioTwo::SD => {
                None
                },
                ScenarioTwo::DD => {
                None
                },
                ScenarioTwo::ICI => {
                None
                },
                ScenarioTwo::ICS => {
                None
                },
                ScenarioTwo::ICD => {
                None
                },
                ScenarioTwo::SCI => {
                None
                },
                ScenarioTwo::SCS => {
                None
                },
                ScenarioTwo::SCD => {
                None
                },
                ScenarioTwo::DCI => {
                None
                },
                ScenarioTwo::DCS => {
                None
                },
                ScenarioTwo::DCD => {
                None
            },*/
            _ => (vec![], 0),
        }
    }
}

pub type Two<'a> = Exist<'a, ScenarioTwo>;

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

        let corrector = Two::new(&set, 2);

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

        let corrector = Two::new(&set, 2);

        assert_eq!(refe, corrector.correct(read).as_slice());
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

        let corrector = Two::new(&set, 2);

        assert_eq!(refe, corrector.correct(read).as_slice());
        assert_eq!(refe, corrector.correct(refe).as_slice());
    }

    #[test]
    #[ignore]
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

        let corrector = Two::new(&set, 2);

        assert_eq!(refe, corrector.correct(read).as_slice());
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

        let corrector = Two::new(&set, 2);

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

        let corrector = Two::new(&set, 2);

        assert_eq!(refe, corrector.correct(read).as_slice());
        assert_eq!(refe, corrector.correct(refe).as_slice());
    }

    #[test]
    #[ignore]
    fn ciic() {
        init();

        let refe = b"GATACATGGACACTAGTATG";
        //           ||||||||||
        let read = b"GATACATGGATTCACTAGTATG";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = Two::new(&set, 2);

        assert_eq!(refe, corrector.correct(read).as_slice()); // test correction work
        assert_eq!(refe, corrector.correct(refe).as_slice()); // test not overcorrection
    }

    #[test]
    fn cisc() {
        init();

        let refe = b"GATACATGGACACTAGTATG";
        //           |||||||||| \\\\\\\\\\
        let read = b"GATACATGGATGACTAGTATG";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = Two::new(&set, 2);

        assert_eq!(refe, corrector.correct(read).as_slice()); // test correction work
        assert_eq!(refe, corrector.correct(refe).as_slice()); // test not overcorrection
    }
}
