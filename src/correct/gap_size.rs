/*
Copyright (c) 2020 Pierre Marijon <pmarijon@hhu.de>

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

pub struct GapSize<'a> {
    valid_kmer: &'a pcon::solid::Solid,
    graph: graph::Graph<'a>,
    one: one::One<'a>,
}

impl<'a> GapSize<'a> {
    pub fn new(valid_kmer: &'a pcon::solid::Solid, c: u8) -> Self {
        Self {
            valid_kmer,
            graph: graph::Graph::new(valid_kmer),
            one: one::One::new(valid_kmer, c),
        }
    }

    pub fn ins_sub_correction(&mut self, kmer: u64, gap_size: usize) -> Option<(Vec<u8>, usize)> {
        let mut alts = alt_nucs(self.valid_kmer, kmer);

        if alts.len() != 1 {
            debug!("not one alts {:?}", alts);
            return None;
        }

        let mut corr = add_nuc_to_end(kmer >> 2, alts[0], self.k());
        let mut local_corr = vec![cocktail::kmer::bit2nuc(alts[0])];

        for i in 0..gap_size {
            debug!("kmer {:?}", cocktail::kmer::kmer2seq(corr, self.k()));

            alts = next_nucs(self.valid_kmer, corr);

            if alts.len() != 1 {
                debug!("failled multiple successor {:?} i: {}", alts, i);
                return None;
            }

            corr = add_nuc_to_end(corr, alts[0], self.k());

            local_corr.push(cocktail::kmer::bit2nuc(alts[0]));
        }

        let offset = local_corr.len();
        Some((local_corr, offset))
    }
}

impl<'a> Corrector for GapSize<'a> {
    fn valid_kmer(&self) -> &pcon::solid::Solid {
        self.valid_kmer
    }

    fn correct_error(&mut self, kmer: u64, seq: &[u8]) -> Option<(Vec<u8>, usize)> {
        let (error_len, _first_correct_kmer) = error_len(&seq, kmer, self.valid_kmer());

        debug!("error_len {}", error_len);
        if error_len < self.k() as usize {
            self.graph.correct_error(kmer, seq) // we can avoid a second compute of error_len
        } else if error_len == self.k() as usize {
            self.one.correct_error(kmer, seq)
        } else {
            self.ins_sub_correction(kmer, error_len - self.k() as usize)
        }
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
    fn csc() {
        init();

        let refe = b"AGCGTATCTT";
        //           ||||| |||||
        let read = b"AGCGTTTCTT";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        let mut corrector = GapSize::new(&data, 2);

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

        let mut corrector = GapSize::new(&data, 2);

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

        let mut corrector = GapSize::new(&data, 2);

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

        let mut corrector = GapSize::new(&data, 2);

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

        let mut corrector = GapSize::new(&data, 2);

        assert_eq!(refe, corrector.correct(read).as_slice());
        assert_eq!(refe, corrector.correct(refe).as_slice());
    }

    #[test]
    fn cddc() {
        init();

        let refe = b"CAAAGCATTTTT";
        //           |||||
        let read = b"CAAAGTTTTT";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        let mut corrector = GapSize::new(&data, 2);

        assert_eq!(refe, corrector.correct(read).as_slice());
        assert_eq!(refe, corrector.correct(refe).as_slice());
    }

    #[test]
    fn cic() {
        init();

        let refe = b"GGATAACTCT";
        //           |||||
        let read = b"GGATATACTCT";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        let mut corrector = GapSize::new(&data, 2);

        assert_eq!(refe, corrector.correct(read).as_slice()); // test correction work
        assert_eq!(refe, corrector.correct(refe).as_slice()); // test not overcorrection
    }

    #[test]
    fn ciic() {
        init();

        let refe = b"GCGTAAATGGAT";
        //           ||||||
        let read = b"GCGTAATTATGGAT";

        let mut data: pcon::solid::Solid = pcon::solid::Solid::new(5);

        for kmer in cocktail::tokenizer::Tokenizer::new(refe, 5) {
            data.set(kmer, true);
        }

        let mut corrector = GapSize::new(&data, 2);

        assert_eq!(refe, corrector.correct(read).as_slice()); // test correction work
        assert_eq!(refe, corrector.correct(refe).as_slice()); // test not overcorrection
    }
}
