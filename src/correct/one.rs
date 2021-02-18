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

#[derive(Debug)]
enum Scenario {
    Deletion,
    Substitution,
    Insertion,
}

impl Scenario {
    fn check(&self) -> u8 {
        match self {
            Scenario::Deletion => 0,
            Scenario::Substitution => 1,
            Scenario::Insertion => 2,
        }
    }

    fn post(&self) -> usize {
        match self {
            Scenario::Deletion => 0,
            Scenario::Substitution => 1,
            Scenario::Insertion => 1,
        }
    }

    fn name(&self) -> &str {
        match self {
            Scenario::Deletion => "Deletion",
            Scenario::Substitution => "Substitution",
            Scenario::Insertion => "Insertion",
        }
    }
}

pub struct One<'a> {
    valid_kmer: &'a set::BoxKmerSet<'a>,
    c: u8,
}

impl<'a> One<'a> {
    pub fn new(valid_kmer: &'a set::BoxKmerSet<'a>, c: u8) -> Self {
        Self { valid_kmer, c }
    }

    fn evaluate_scenario(
        &self,
        scenari: Scenario,
        corr: u64,
        seq: &[u8],
        scenario: &mut Vec<(Scenario, u64, Vec<u8>)>,
        alts: Vec<u8>,
    ) {
        if get_score(
            self.valid_kmer,
            corr,
            &seq[scenari.check() as usize..(self.c + scenari.check()) as usize],
        ) == self.c as usize
        {
            debug!("it's a {}", scenari.name());
            scenario.push((scenari, corr, alts));
        }
    }

    fn get_scenario(&self, kmer: u64, seq: &[u8]) -> Vec<(Scenario, u64, Vec<u8>)> {
        let mut scenario: Vec<(Scenario, u64, Vec<u8>)> = Vec::new();

        let alts = alt_nucs(self.valid_kmer, kmer);

        if alts.len() != 1 {
            debug!("not one alts {:?}", alts);
            return scenario;
        }
        debug!("one alts {:?}", alts);

        let corr = add_nuc_to_end(kmer >> 2, alts[0], self.k());

        if (self.c + Scenario::Insertion.check()) as usize <= seq.len() {
            debug!("Test deletion");
            self.evaluate_scenario(
                Scenario::Deletion,
                corr,
                seq,
                &mut scenario,
                vec![cocktail::kmer::bit2nuc(alts[0])],
            );
            debug!("Test substitution");
            self.evaluate_scenario(
                Scenario::Substitution,
                corr,
                seq,
                &mut scenario,
                vec![cocktail::kmer::bit2nuc(alts[0])],
            );
            debug!("Test insertion");
            self.evaluate_scenario(Scenario::Insertion, corr, seq, &mut scenario, vec![]);
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
            Some((scenario[0].2.clone(), scenario[0].0.post()))
        } else if (self.c + Scenario::Insertion.check() + 1) as usize <= seq.len() {
            scenario.retain(|x| {
                last_is_valid(
                    self.valid_kmer,
                    x.1,
                    &seq[(x.0.check() as usize)..(self.c + x.0.check() + 1) as usize],
                )
            });

            if scenario.len() == 1 {
                debug!("multi scenario one is better {}", scenario[0].0.name());
                Some((scenario[0].2.clone(), scenario[0].0.post()))
            } else {
                debug!("multi scenario no better");
                None
            }
        } else {
            None
        }
    }
}

fn get_score(valid_kmer: &set::BoxKmerSet, mut kmer: u64, nucs: &[u8]) -> usize {
    let mut score = 0;

    for nuc in nucs {
        debug!("before {}", cocktail::kmer::kmer2seq(kmer, valid_kmer.k()));
        kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(*nuc), valid_kmer.k());
        debug!(
            "after  {} {}",
            cocktail::kmer::kmer2seq(kmer, valid_kmer.k()),
            valid_kmer.get(kmer)
        );
        if valid_kmer.get(kmer) {
            score += 1
        } else {
            break;
        }
    }

    score
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
