//! Pcon set

/* std use */

/* crates use */
#[cfg(feature = "parallel")]
use rayon::prelude::*;

/* project use */
use crate::error;
use crate::set;

pub struct Pcon {
    set: pcon::solid::Solid,
}

impl Pcon {
    pub fn from_pcon_solid<R>(input: R) -> error::Result<Self>
    where
        R: std::io::BufRead,
    {
        Ok(Self {
            set: pcon::solid::Solid::from_stream(input)?,
        })
    }

    pub fn from_csv<R>(input: R, k: u8) -> error::Result<Self>
    where
        R: std::io::BufRead,
    {
        let mut set = pcon::solid::Solid::new(k);

        let mut reader = csv::Reader::from_reader(input);
        for result in reader.byte_records() {
            let record = result?;

            set.set(
                cocktail::kmer::seq2bit(record.get(0).ok_or(error::Error::CsvMissingFirstColumn)?),
                true,
            );
        }

        Ok(Self { set })
    }

    #[cfg(not(feature = "parallel"))]
    pub fn from_fasta<R>(input: R, k: u8) -> Self
    where
        R: std::io::BufRead,
    {
        let mut set = pcon::solid::Solid::new(k);

        let mut reader = noodles::fasta::Reader::new(input);
        let mut records = reader.records();

        while let Some(Ok(record)) = records.next() {
            if record.sequence().len() >= k as usize {
                let kmerizer = cocktail::tokenizer::Canonical::new(record.sequence().as_ref(), k);

                for canonical in kmerizer {
                    set.set(canonical, true);
                }
            }
        }

        Self { set }
    }

    #[cfg(feature = "parallel")]
    pub fn from_fasta<R>(input: R, k: u8) -> Self
    where
        R: std::io::BufRead,
    {
        let mut set = pcon::solid::Solid::new(k);

        let mut reader = noodles::fasta::Reader::new(input);

        let mut iter = reader.records();
        let mut records = Vec::with_capacity(8192);

        let mut end = true;
        while end {
            log::info!("Start populate buffer");
            end = crate::populate_buffer(&mut iter, &mut records, 8192);
            log::info!("End populate buffer {}", records.len());

            set.extend(
                records
                    .par_iter()
                    .filter(|record| record.sequence().len() >= k as usize)
                    .map(|record| {
                        let mut set = pcon::solid::Solid::new(k);
                        let kmerizer =
                            cocktail::tokenizer::Canonical::new(record.sequence().as_ref(), k);
                        for canonical in kmerizer {
                            set.set(canonical, true);
                        }
                        set
                    })
                    .reduce(
                        || pcon::solid::Solid::new(k),
                        |mut x, y| {
                            x.extend(y);
                            x
                        },
                    ),
            );
        }

        Self { set }
    }

    #[cfg(not(feature = "parallel"))]
    pub fn from_fastq<R>(input: R, k: u8) -> Self
    where
        R: std::io::BufRead,
    {
        let mut set = pcon::solid::Solid::new(k);

        let mut reader = noodles::fastq::Reader::new(input);
        let mut records = reader.records();

        while let Some(Ok(record)) = records.next() {
            if record.sequence().len() >= k as usize {
                let kmerizer = cocktail::tokenizer::Canonical::new(record.sequence(), k);

                for canonical in kmerizer {
                    set.set(canonical, true);
                }
            }
        }

        Self { set }
    }

    #[cfg(feature = "parallel")]
    pub fn from_fastq<R>(input: R, k: u8) -> Self
    where
        R: std::io::BufRead,
    {
        let mut set = pcon::solid::Solid::new(k);

        let mut reader = noodles::fastq::Reader::new(input);

        let mut iter = reader.records();
        let mut records = Vec::with_capacity(8192);

        let mut end = true;
        while end {
            log::info!("Start populate buffer");
            end = crate::populate_bufferq(&mut iter, &mut records, 8192);
            log::info!("End populate buffer {}", records.len());

            set.extend(
                records
                    .par_iter()
                    .filter(|record| record.sequence().len() >= k as usize)
                    .map(|record| {
                        let mut set = pcon::solid::Solid::new(k);
                        let kmerizer =
                            cocktail::tokenizer::Canonical::new(record.sequence().as_ref(), k);
                        for canonical in kmerizer {
                            set.set(canonical, true);
                        }
                        set
                    })
                    .reduce(
                        || pcon::solid::Solid::new(k),
                        |mut x, y| {
                            x.extend(y);
                            x
                        },
                    ),
            );
        }

        Self { set }
    }

    pub fn new(set: pcon::solid::Solid) -> Self {
        Self { set }
    }
}

impl set::KmerSet for Pcon {
    fn get(&self, kmer: u64) -> bool {
        self.set.get(kmer)
    }

    fn k(&self) -> u8 {
        self.set.k()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    static SEQ: &[u8] = b"ACGTGGGAATTGTGGCCACATCACGAGGTCCTGCGTATTGACGACTGTAAAGCGAGTGGCCGTGGAATTTCAAGCTCAATTAGCCGAACCAATCCGCCTA";

    #[test]
    fn canonical() {
        let mut solid = pcon::solid::Solid::new(11);
        for cano in cocktail::tokenizer::Canonical::new(SEQ, 11) {
            solid.set(cano, true);
        }

        let set: crate::set::BoxKmerSet = Box::new(Pcon::new(solid));

        for cano in cocktail::tokenizer::Canonical::new(SEQ, 11) {
            assert!(set.get(cano))
        }
    }

    #[test]
    fn forward() {
        let mut solid = pcon::solid::Solid::new(11);
        for cano in cocktail::tokenizer::Canonical::new(SEQ, 11) {
            solid.set(cano, true);
        }

        let set: crate::set::BoxKmerSet = Box::new(Pcon::new(solid));

        for kmer in cocktail::tokenizer::Tokenizer::new(SEQ, 11) {
            assert!(set.get(kmer))
        }
    }

    #[test]
    fn absence() {
        let mut solid = pcon::solid::Solid::new(11);
        for cano in cocktail::tokenizer::Canonical::new(SEQ, 11) {
            solid.set(cano, true);
        }

        let set: crate::set::BoxKmerSet = Box::new(Pcon::new(solid));

        assert!(!set.get(0));
    }

    #[test]
    fn k() {
        let mut solid = pcon::solid::Solid::new(11);
        for cano in cocktail::tokenizer::Canonical::new(SEQ, 11) {
            solid.set(cano, true);
        }

        let set: crate::set::BoxKmerSet = Box::new(Pcon::new(solid));

        assert_eq!(set.k(), 11);
    }
}
