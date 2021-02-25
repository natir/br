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
    DCD(usize, u8),
}

impl Scenario for ScenarioTwo {
    fn init(&self, c: usize, k: u8) -> Self {
        match self {
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
            ScenarioTwo::DCD(_, _) => ScenarioTwo::DCD(c, k),
        }
    }

    fn c(&self) -> usize {
        match self {
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
            ScenarioTwo::II(_, _) => Some((kmer, 3)), // kmer not change check from base 3
            ScenarioTwo::IS(_, _) => Some((kmer, 2)), // kmer not change check from base 2
            ScenarioTwo::SS(_, k) => {
                if seq.is_empty() {
                    None
                } else {
                    kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(seq[1]), *k);
                    if valid_kmer.get(kmer) {
                        None
                    } else {
                        let alts = alt_nucs(valid_kmer, kmer);
                        if alts.len() != 1 {
                            None
                        } else {
                            Some((add_nuc_to_end(kmer >> 2, alts[0], *k), 2))
                        }
                    }
                }
            }
            ScenarioTwo::SD(_, k) => {
                if seq.is_empty() {
                    None
                } else {
                    let alts = alt_nucs(valid_kmer, kmer << 2);
                    if alts.len() != 1 {
                        None
                    } else {
                        Some((add_nuc_to_end(kmer, alts[0], *k), 1))
                    }
                }
            }
            ScenarioTwo::DD(_, k) => {
                let alts = alt_nucs(valid_kmer, kmer << 2);
                if alts.len() != 1 {
                    None
                } else {
                    Some((add_nuc_to_end(kmer, alts[0], *k), 0))
                }
            }
            ScenarioTwo::ICI(_, k) => {
                // If seq is to short we can't apply
                if seq.len() < 4 {
                    None
                } else {
                    let corr = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(seq[3]), *k); // If kmer with thrid base is valid is an ICI

                    if valid_kmer.get(corr) {
                        Some((corr, 4))
                    } else {
                        None
                    }
                }
            }
            ScenarioTwo::ICS(_, k) => {
                // If seq is to short we can't apply
                if seq.len() < 4 {
                    None
                } else {
                    kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(seq[1]), *k);
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
            ScenarioTwo::ICD(_, k) => {
                // If seq is to short we can't apply
                if seq.len() < 4 {
                    None
                } else {
                    let second = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(seq[2]), *k);

                    let alts = alt_nucs(valid_kmer, second << 2);
                    if alts.len() != 1 {
                        None
                    } else {
                        Some((add_nuc_to_end(second, alts[0], *k), 3))
                    }
                }
            }
            ScenarioTwo::SCI(_, k) => {
                if seq.len() < 4 {
                    None
                } else {
                    kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(seq[1]), *k);
                    kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(seq[3]), *k);

                    Some((kmer, 4))
                }
            }
            ScenarioTwo::SCS(_, k) => {
                if seq.len() < 3 {
                    None
                } else {
                    kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(seq[1]), *k);
                    if valid_kmer.get(kmer) {
                        kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(seq[2]), *k);

                        if !valid_kmer.get(kmer) {
                            let alts = alt_nucs(valid_kmer, kmer);

                            if alts.len() == 1 {
                                Some((add_nuc_to_end(kmer >> 2, alts[0], *k), 3))
                            } else {
                                None
                            }
                        } else {
                            None
                        }
                    } else {
                        None
                    }
                }
            }
            ScenarioTwo::SCD(_, k) => {
                if seq.len() < 2 {
                    None
                } else {
                    kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(seq[1]), *k);

                    let alts = alt_nucs(valid_kmer, kmer << 2);

                    if alts.len() != 1 {
                        None
                    } else {
                        Some((add_nuc_to_end(kmer, alts[0], *k), 2))
                    }
                }
            }
            ScenarioTwo::DCI(_, k) => {
                if seq.len() < 4 {
                    None
                } else {
                    kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(seq[1]), *k);
                    kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(seq[3]), *k);

                    Some((kmer, 4))
                }
            }
            ScenarioTwo::DCD(_, k) => {
                if seq.len() < 2 {
                    None
                } else {
                    kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(seq[0]), *k);

                    let alts = alt_nucs(valid_kmer, kmer << 2);
                    if alts.len() != 1 {
                        None
                    } else {
                        Some((add_nuc_to_end(kmer, alts[0], *k), 1))
                    }
                }
            }
        }
    }

    fn correct(&self, valid_kmer: &set::BoxKmerSet, kmer: u64, seq: &[u8]) -> (Vec<u8>, usize) {
        match self {
            ScenarioTwo::II(_, _) => (vec![cocktail::kmer::bit2nuc(kmer & 0b11)], 2),
            ScenarioTwo::IS(_, _) => (vec![cocktail::kmer::bit2nuc(kmer & 0b11)], 2),
            ScenarioTwo::SS(_, _) | ScenarioTwo::SD(_, _) | ScenarioTwo::DD(_, _) => {
                let (corr, offset) = self
                    .apply(valid_kmer, kmer, seq)
                    .expect("we can't failled her");

                (
                    vec![
                        cocktail::kmer::bit2nuc((corr & 0b1100) >> 2),
                        cocktail::kmer::bit2nuc(corr & 0b11),
                    ],
                    offset,
                )
            }
            ScenarioTwo::ICI(_, _) => (vec![cocktail::kmer::bit2nuc(kmer & 0b11)], 3),
            ScenarioTwo::ICD(_, _) => {
                let (corr, offset) = self
                    .apply(valid_kmer, kmer, seq)
                    .expect("we can't failled her");

                (
                    vec![
                        cocktail::kmer::bit2nuc((corr & 0b1100) >> 2),
                        cocktail::kmer::bit2nuc(corr & 0b11),
                    ],
                    offset - 1,
                )
            }
            ScenarioTwo::ICS(_, _) => {
                let (corr, offset) = self
                    .apply(valid_kmer, kmer, seq)
                    .expect("we can't failled her");

                (
                    vec![
                        cocktail::kmer::bit2nuc((corr & 0b1100) >> 2),
                        cocktail::kmer::bit2nuc(corr & 0b11),
                    ],
                    offset + 1,
                )
            }
            ScenarioTwo::SCI(_, _)
            | ScenarioTwo::SCS(_, _)
            | ScenarioTwo::SCD(_, _)
            | ScenarioTwo::DCD(_, _) => {
                let (corr, offset) = self
                    .apply(valid_kmer, kmer, seq)
                    .expect("we can't failled her");

                (
                    vec![
                        cocktail::kmer::bit2nuc((corr & 0b110000) >> 4),
                        cocktail::kmer::bit2nuc((corr & 0b1100) >> 2),
                        cocktail::kmer::bit2nuc(corr & 0b11),
                    ],
                    offset,
                )
            }
            /*            ScenarioTwo::DCD => {
                           None
                       },
            */
            _ => (vec![], 1),
        }
    }
}

pub type Two<'a> = Exist<'a, ScenarioTwo>;

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
    fn ciic() {
        init();

        let refe = filter(b"GATACATGGA--CACTAGTATG");
        //                  ||||||||||  ||||||||||
        let read = filter(b"GATACATGGATTCACTAGTATG");

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(&refe, 5) {
            data.set(kmer, true);
        }

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = Two::new(&set, 2);

        assert_eq!(refe, corrector.correct(&read).as_slice()); // test correction work
        assert_eq!(refe, corrector.correct(&refe).as_slice()); // test not overcorrection
    }

    #[test]
    fn cisc() {
        init();

        let refe = filter(b"GATACATGGA-CACTAGTATG");
        //                  ||||||||||  |||||||||
        let read = filter(b"GATACATGGATGACTAGTATG");

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(7);

        for kmer in cocktail::tokenizer::Tokenizer::new(&refe, 7) {
            data.set(kmer, true);
        }

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = Two::new(&set, 2);

        assert_eq!(refe, corrector.correct(&read).as_slice()); // test correction work
        assert_eq!(refe, corrector.correct(&refe).as_slice()); // test not overcorrection
    }

    #[test]
    fn cssc() {
        init();

        let refe = b"TCGTTATTCGGTGGACTCCT";
        //           ||||||||||  ||||||||
        let read = b"TCGTTATTCGAAGGACTCCT";

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
    fn csdc() {
        init();

        let refe = b"AACAGCTGAATCTACCATTG";
        //           |||||||||| /////////
        let read = b"AACAGCTGAAGTACCATTG";

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
    fn cddc() {
        init();

        let refe = b"TGCCGTAGGCCATTGCGGCT";
        //           ||||||||||  ||||||||
        let read = b"TGCCGTAGGC--TTGCGGCT";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(7);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 7) {
            data.set(kmer, true);
        }

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = Two::new(&set, 2);

        assert_eq!(refe, corrector.correct(&filter(read)).as_slice()); // test correction work
        assert_eq!(refe, corrector.correct(refe).as_slice()); // test not overcorrection
    }

    #[test]
    fn cicic() {
        init();

        let refe = filter(b"ATAGTAACGG-A-CACACTT");
        //                  |||||||||| | |||||||
        let read = filter(b"ATAGTAACGGAAGCACACTT");

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(7);

        for kmer in cocktail::tokenizer::Tokenizer::new(&refe, 7) {
            data.set(kmer, true);
        }

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = Two::new(&set, 3);

        assert_eq!(&refe, corrector.correct(&read).as_slice()); // test correction work
        assert_eq!(&refe, corrector.correct(&refe).as_slice()); // test not overcorrection
    }

    #[test]
    fn cicsc() {
        init();

        let refe = filter(b"GAGCCCAGAG-CGATATTCT");
        //                  |||||||||| | |||||||
        let read = filter(b"GAGCCCAGAGACTATATTCT");

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(7);

        for kmer in cocktail::tokenizer::Tokenizer::new(&refe, 7) {
            data.set(kmer, true);
        }

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = Two::new(&set, 2);

        assert_eq!(&refe, corrector.correct(&read).as_slice()); // test correction work
        assert_eq!(&refe, corrector.correct(&refe).as_slice()); // test not overcorrection
    }

    #[test]
    fn cicdc() {
        init();

        let refe = filter(b"TCGAAAGCAT-GGGTACGTT");
        //                  |||||||||| | |||||||
        let read = filter(b"TCGAAAGCATAG-GTACGTT");

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(7);

        for kmer in cocktail::tokenizer::Tokenizer::new(&refe, 7) {
            data.set(kmer, true);
        }

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = Two::new(&set, 2);

        assert_eq!(&refe, corrector.correct(&read).as_slice()); // test correction work
        assert_eq!(&refe, corrector.correct(&refe).as_slice()); // test not overcorrection
    }

    #[test]
    fn cscic() {
        init();

        let refe = filter(b"AAGGATGCATCG-ACTCAAG");
        //                  |||||||||| | |||||||
        let read = filter(b"AAGGATGCATGGAACTCAAG");

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(7);

        for kmer in cocktail::tokenizer::Tokenizer::new(&refe, 7) {
            data.set(kmer, true);
        }

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = Two::new(&set, 2);

        assert_eq!(&refe, corrector.correct(&read).as_slice()); // test correction work
        assert_eq!(&refe, corrector.correct(&refe).as_slice()); // test not overcorrection
    }

    #[test]
    fn cscsc() {
        init();

        let refe = filter(b"ACACGTGCGCTTGGAGGTAC");
        //                  |||||||||| | |||||||
        let read = filter(b"ACACGTGCGCATCGAGGTAC");

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(7);

        for kmer in cocktail::tokenizer::Tokenizer::new(&refe, 7) {
            data.set(kmer, true);
        }

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = Two::new(&set, 2);

        assert_eq!(&refe, corrector.correct(&read).as_slice()); // test correction work
        assert_eq!(&refe, corrector.correct(&refe).as_slice()); // test not overcorrection
    }

    #[test]
    fn cscdc() {
        init();

        let refe = filter(b"TATGCTCTGCGTAATCATAG");
        //                  |||||||||| | |||||||
        let read = filter(b"TATGCTCTGCAT-ATCATAG");

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(7);

        for kmer in cocktail::tokenizer::Tokenizer::new(&refe, 7) {
            data.set(kmer, true);
        }

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = Two::new(&set, 2);

        assert_eq!(&refe, corrector.correct(&read).as_slice()); // test correction work
        assert_eq!(&refe, corrector.correct(&refe).as_slice()); // test not overcorrection
    }

    #[test]
    fn cdcic() {
        init();

        let refe = filter(b"GCTTCGTGATAG-TACGCTT");
        //                  |||||||||| | |||||||
        let read = filter(b"GCTTCGTGAT-GATACGCTT");

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(7);

        for kmer in cocktail::tokenizer::Tokenizer::new(&refe, 7) {
            data.set(kmer, true);
        }

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = Two::new(&set, 2);

        assert_eq!(&refe, corrector.correct(&read).as_slice()); // test correction work
        assert_eq!(&refe, corrector.correct(&refe).as_slice()); // test not overcorrection
    }

    #[test]
    fn cdcsc() {
        init();

        let refe = filter(b"GGACCTGATCACGTCAATTA");
        //                  |||||||||| | |||||||
        let read = filter(b"GGACCTGATC-CCTCAATTA");

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(7);

        for kmer in cocktail::tokenizer::Tokenizer::new(&refe, 7) {
            data.set(kmer, true);
        }

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = Two::new(&set, 2);

        assert_eq!(&refe, corrector.correct(&read).as_slice()); // test correction work
        assert_eq!(&refe, corrector.correct(&refe).as_slice()); // test not overcorrection
    }

    #[test]
    fn cdcdc() {
        init();

        let refe = filter(b"GGAATACGTGCGTTGGGTAA");
        //                  |||||||||| | |||||||
        let read = filter(b"GGAATACGTG-G-TGGGTAA");

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(7);

        for kmer in cocktail::tokenizer::Tokenizer::new(&refe, 7) {
            data.set(kmer, true);
        }

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(data));

        let corrector = Two::new(&set, 2);

        assert_eq!(&refe, corrector.correct(&read).as_slice()); // test correction work
        assert_eq!(&refe, corrector.correct(&refe).as_slice()); // test not overcorrection
    }
}
