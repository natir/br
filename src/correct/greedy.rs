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
use log::{debug, trace};

/* local use */
use crate::correct::*;

struct Score;

impl bio::alignment::pairwise::MatchFunc for Score {
    fn score(&self, a: u8, b: u8) -> i32 {
        if a == b {
            1
        } else {
            -1
        }
    }
}

pub struct Greedy<'a> {
    valid_kmer: &'a pcon::solid::Solid,
    max_search: u8,
    nb_validate: u8,
}

impl<'a> Greedy<'a> {
    pub fn new(valid_kmer: &'a pcon::solid::Solid, max_search: u8, nb_validate: u8) -> Self {
        Self {
            valid_kmer,
            max_search,
            nb_validate,
        }
    }

    fn match_alignement(&self, before_seq: Vec<u8>, read: &[u8], corr: &[u8]) -> Option<i64> {
        let mut r = before_seq.clone();
        r.extend_from_slice(read);

        let mut c = before_seq.clone();
        c.extend_from_slice(corr);

        let mut aligner =
            bio::alignment::pairwise::Aligner::with_capacity(10, 10, -1, -1, Score {});
        let alignment = aligner.global(r.as_slice(), c.as_slice());

        trace!("{}", alignment.pretty(r.as_slice(), c.as_slice()));

        let mut offset = 0;
        for ops in alignment.operations[before_seq.len()..].windows(2) {
            match ops[0] {
                bio::alignment::AlignmentOperation::Del => offset -= 1,
                bio::alignment::AlignmentOperation::Ins => offset += 1,
                _ => (),
            }

            if ops[0] == bio::alignment::AlignmentOperation::Match && ops[0] == ops[1] {
                return Some(offset);
            }
        }

        None
    }

    fn follow_graph(&self, mut kmer: u64) -> Option<(u8, u64)> {
        let alts = next_nucs(self.valid_kmer(), kmer);

        if alts.len() != 1 {
            debug!("failled branching node {:?}", alts);
            return None;
        }

        kmer = add_nuc_to_end(kmer, alts[0], self.k());

        Some((cocktail::kmer::bit2nuc(alts[0]), kmer))
    }

    fn check_next_kmers(&self, mut kmer: u64, seq: &[u8]) -> bool {
        if seq.len() < self.nb_validate as usize {
            return false;
        }

        for nuc in &seq[..self.nb_validate as usize] {
            kmer = add_nuc_to_end(kmer, cocktail::kmer::nuc2bit(*nuc), self.k());
            if !self.valid_kmer.get(kmer) {
                return false;
            }
        }

        true
    }
}

impl<'a> Corrector for Greedy<'a> {
    fn k(&self) -> u8 {
        self.valid_kmer.k
    }

    fn valid_kmer(&self) -> &pcon::solid::Solid {
        self.valid_kmer
    }

    fn correct_error(&self, mut kmer: u64, seq: &[u8]) -> Option<(Vec<u8>, usize)> {
        let alts = alt_nucs(self.valid_kmer(), kmer);
        if alts.len() != 1 {
            debug!("failled multiple successor {:?}", alts);
            return None;
        }

        let mut local_corr = Vec::new();
        let before_seq = cocktail::kmer::kmer2seq(kmer >> 2, self.k() - 1)
            .as_bytes()
            .to_vec();

        kmer = add_nuc_to_end(kmer >> 2, alts[0], self.k());

        local_corr.push(cocktail::kmer::bit2nuc(alts[0]));

        for i in 0..(self.max_search as usize) {
            if let Some((base, new_kmer)) = self.follow_graph(kmer) {
                local_corr.push(base);
                kmer = new_kmer;
            }

            if seq.len() < i as usize {
                return None;
            }

            if let Some(off) = self.match_alignement(before_seq.clone(), &seq[..i], &local_corr) {
                if self.check_next_kmers(kmer, &seq[i..]) {
                    let offset: usize = (local_corr.len() as i64 + off) as usize;
                    return Some((local_corr, offset));
                }
            }
        }

        None
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    fn init() {
        let _ = env_logger::builder()
            .is_test(true)
            .filter_level(log::LevelFilter::Trace)
            .try_init();
    }

    #[test]
    fn branching_path_csc() {
        init();

        let refe = b"TCTTTATTTTC";
        //           ||||| |||||
        let read = b"TCTTTGTTTTC";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        data.set(cocktail::kmer::seq2bit(b"CTTTT"), true);

        let corrector = Greedy::new(&data, 7, 2);

        assert_eq!(read, corrector.correct(read).as_slice()); // test correction work
        assert_eq!(refe, corrector.correct(refe).as_slice()); // test not overcorrection
    }

    #[test]
    fn branching_path_cdc() {
        init();

        let refe = b"GATACATGGACACTAGTATG";
        //           ||||||||||
        let read = b"GATACATGGAACTAGTATG";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        data.set(cocktail::kmer::seq2bit(b"GGACT"), true);

        let corrector = Greedy::new(&data, 7, 2);

        assert_eq!(read, corrector.correct(read).as_slice());
        assert_eq!(refe, corrector.correct(refe).as_slice());
    }

    #[test]
    fn branching_path_cic() {
        init();

        let refe = b"GATACATGGACACTAGTATG";
        //           ||||||||||
        let read = b"GATACATGGATCACTAGTATG";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 7) {
            data.set(kmer, true);
        }

        data.set(cocktail::kmer::seq2bit(b"GGACT"), true);

        let corrector = Greedy::new(&data, 7, 2);

        assert_eq!(read, corrector.correct(read).as_slice()); // test correction work
        assert_eq!(refe, corrector.correct(refe).as_slice()); // test not overcorrection
    }

    #[test]
    fn csc() {
        init();

        let refe = b"TCTTTATTTTC";
        //           ||||| |||||
        let read = b"TCTTTGTTTTC";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        let corrector = Greedy::new(&data, 7, 2);

        assert_eq!(refe, corrector.correct(read).as_slice()); // test correction work
        assert_eq!(refe, corrector.correct(refe).as_slice()); // test not overcorrection
    }

    #[test]
    fn cssc() {
        init();

        let refe = b"TCTCTAATCTTC";
        //           |||||  |||||
        let read = b"TCTCTGGTCTTC";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        let corrector = Greedy::new(&data, 7, 2);

        assert_eq!(refe, corrector.correct(read).as_slice()); // test correction work
        assert_eq!(refe, corrector.correct(refe).as_slice()); // test not overcorrection
    }

    #[test]
    fn csssc() {
        init();

        let refe = b"TCTCTAAATCTTC";
        //           |||||  |||||
        let read = b"TCTCTGGGTCTTC";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        let corrector = Greedy::new(&data, 7, 2);

        assert_eq!(refe, corrector.correct(read).as_slice()); // test correction work
        assert_eq!(refe, corrector.correct(refe).as_slice()); // test not overcorrect
    }

    #[test]
    fn cscsc() {
        init();

        let refe = b"TCTTTACATTTTT";
        //           ||||| | |||||
        let read = b"TCTTTGCGTTTTT";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        let corrector = Greedy::new(&data, 7, 2);

        assert_eq!(refe, corrector.correct(read).as_slice()); // test correction work
        assert_eq!(refe, corrector.correct(refe).as_slice()); // test not overcorrection
    }

    #[test]
    fn cdc() {
        init();

        let refe = b"GATACATGGACACTAGTATG";
        //           ||||||||||
        let read = b"GATACATGGAACTAGTATG";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        let corrector = Greedy::new(&data, 7, 2);

        assert_eq!(refe, corrector.correct(read).as_slice());
        assert_eq!(refe, corrector.correct(refe).as_slice());
    }

    #[test]
    fn cddc() {
        init();

        let refe = b"CAAAGCATTTTTT";
        //           |||||
        let read = b"CAAAGTTTTTT";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        let corrector = Greedy::new(&data, 7, 2);

        assert_eq!(refe, corrector.correct(read).as_slice());
        assert_eq!(refe, corrector.correct(refe).as_slice());
    }

    #[test]
    fn cic() {
        init();

        let refe = b"GATACATGGACACTAGTATG";
        //           ||||||||||
        let read = b"GATACATGGATCACTAGTATG";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        let corrector = Greedy::new(&data, 7, 2);

        assert_eq!(refe, corrector.correct(read).as_slice()); // test correction work
        assert_eq!(refe, corrector.correct(refe).as_slice()); // test not overcorrection
    }

    #[test]
    fn ciic() {
        init();

        let refe = b"GATACATGGACACTAGTATG";
        //           ||||||||||
        let read = b"GATACATGGATTCACTAGTATG";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        let corrector = Greedy::new(&data, 7, 2);

        assert_eq!(refe, corrector.correct(read).as_slice()); // test correction work
        assert_eq!(refe, corrector.correct(refe).as_slice()); // test not overcorrection
    }
}
