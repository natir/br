/* crate use */
use log::debug;
use strum::IntoEnumIterator;

/* crate use */
use crate::correct::*;
use crate::set;

///////////////////////////////////////////
// Generic Trait for correction scenario //
///////////////////////////////////////////
pub trait Scenario: std::fmt::Debug + Copy {
    fn init(&self, c: usize, k: u8) -> Self;

    fn c(&self) -> usize;

    fn apply(&self, valid_kmer: &set::BoxKmerSet, kmer: u64, seq: &[u8]) -> Option<(u64, usize)>;

    fn correct(&self, kmer: u64, _seq: &[u8]) -> (Vec<u8>, usize);

    fn get_score(&self, valid_kmer: &set::BoxKmerSet, ori: u64, seq: &[u8]) -> usize {
        if let Some((mut kmer, offset)) = self.apply(valid_kmer, ori, seq) {
            if !valid_kmer.get(kmer) {
                return 0;
            }

            if offset + self.c() > seq.len() {
                return 0;
            }

            let mut score = 0;

            for nuc in &seq[offset..offset + self.c()] {
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

    fn one_more(&self, valid_kmer: &set::BoxKmerSet, mut kmer: u64, seq: &[u8]) -> bool {
        let (correct, offset) = self.correct(kmer, seq);
        println!("{:?} correct {:?}, offset {}", self, correct, offset);

        if seq.len() > self.c() + offset + 1 {
            println!("seq[offset.. ] {:?}", &seq[offset..self.c() + offset + 1]);
            println!(
                "correction {}",
                cocktail::kmer::kmer2seq(kmer, valid_kmer.k())
            );
            seq[offset..self.c() + offset + 1].iter().for_each(|nuc| {
                kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(*nuc), valid_kmer.k())
            });
            println!("check {}", cocktail::kmer::kmer2seq(kmer, valid_kmer.k()));

            kmer = add_nuc_to_end(
                kmer,
                cocktail::kmer::nuc2bit(seq[self.c() + offset + 1]),
                valid_kmer.k(),
            );
            valid_kmer.get(kmer)
        } else {
            false
        }
    }
}

//////////////////////////////////////////
// Exsist use scenario to correct error //
//////////////////////////////////////////
pub struct Exist<'a, S>
where
    S: Scenario + IntoEnumIterator,
{
    valid_kmer: &'a set::BoxKmerSet<'a>,
    c: u8,
    _phantom: std::marker::PhantomData<&'a S>,
}

impl<'a, S> Exist<'a, S>
where
    S: Scenario + IntoEnumIterator,
{
    pub fn new(valid_kmer: &'a set::BoxKmerSet, c: u8) -> Self {
        Self {
            valid_kmer,
            c,
            _phantom: std::marker::PhantomData,
        }
    }

    fn get_scenarii(&self, kmer: u64, seq: &[u8]) -> Vec<S> {
        let mut scenarii: Vec<S> = Vec::new();

        for mut scenario in S::iter() {
            scenario = scenario.init(self.c as usize, self.valid_kmer.k());

            if scenario.get_score(self.valid_kmer, kmer, seq) == self.c as usize {
                scenarii.push(scenario)
            }
        }

        scenarii
    }
}

impl<'a, S> Corrector for Exist<'a, S>
where
    S: Scenario + IntoEnumIterator,
{
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
            debug!("one {:?}", scenarii);
            Some(scenarii[0].correct(corr, seq))
        } else {
            debug!("multiple {:?}", scenarii);
            scenarii.retain(|x| x.one_more(self.valid_kmer, corr, seq));
            debug!("multiple {:?}", scenarii);

            if scenarii.len() == 1 {
                Some(scenarii[0].correct(corr, seq))
            } else {
                None
            }
        }
    }
}

pub mod one;
pub mod two;
