//! HashSet

/* std use */

/* crates use */
#[cfg(feature = "parallel")]
use rayon::prelude::*;

/* project use */
use crate::error;
use crate::set;

pub struct Hash {
    set: rustc_hash::FxHashSet<u64>,
    k: u8,
}

impl Hash {
    pub fn from_csv<R>(input: R, k: u8) -> error::Result<Self>
    where
        R: std::io::BufRead,
    {
        let mut set = rustc_hash::FxHashSet::default();

        let mut reader = csv::Reader::from_reader(input);
        for result in reader.byte_records() {
            let record = result?;

            set.insert(cocktail::kmer::canonical(
                cocktail::kmer::seq2bit(record.get(0).ok_or(error::Error::CsvMissingFirstColumn)?),
                k,
            ));
        }

        Ok(Self { set, k })
    }

    #[cfg(not(feature = "parallel"))]
    pub fn from_fasta<R>(input: R, k: u8) -> Self
    where
        R: std::io::BufRead,
    {
        let mut set = rustc_hash::FxHashSet::default();

        let mut reader = noodles::fasta::Reader::new(input);
        let mut records = reader.records();

        while let Some(Ok(record)) = records.next() {
            if record.sequence().len() >= k as usize {
                let kmerizer = cocktail::tokenizer::Canonical::new(record.sequence().as_ref(), k);

                for canonical in kmerizer {
                    set.insert(canonical);
                }
            }
        }

        Self { set, k }
    }

    #[cfg(feature = "parallel")]
    pub fn from_fasta<R>(input: R, k: u8) -> Self
    where
        R: std::io::BufRead,
    {
        let mut set = rustc_hash::FxHashSet::default();

        let mut reader = noodles::fasta::Reader::new(input);
        let mut iter = reader.records();
        let mut records = Vec::with_capacity(8192);

        let mut end = true;
        while end {
            log::info!("Start populate buffer");
            end = populate_buffer(&mut iter, &mut records, 8192);
            log::info!("End populate buffer {}", records.len());

            set.extend(
                records
                    .par_iter()
                    .filter(|record| record.sequence().len() >= k as usize)
                    .map(|record| {
                        let mut set = rustc_hash::FxHashSet::default();
                        let kmerizer =
                            cocktail::tokenizer::Canonical::new(record.sequence().as_ref(), k);
                        for canonical in kmerizer {
                            set.insert(canonical);
                        }
                        set
                    })
                    .reduce(
                        || rustc_hash::FxHashSet::default(),
                        |mut x, y| {
                            x.extend(y.iter());
                            x
                        },
                    ),
            );
        }

        Self { set, k }
    }

    #[cfg(not(feature = "parallel"))]
    pub fn from_fastq<R>(input: R, k: u8) -> Self
    where
        R: std::io::BufRead,
    {
        let mut set = rustc_hash::FxHashSet::default();

        let mut reader = noodles::fasta::Reader::new(input);
        let mut records = reader.records();

        while let Some(Ok(record)) = records.next() {
            if record.sequence().len() >= k as usize {
                let kmerizer = cocktail::tokenizer::Canonical::new(record.sequence().as_ref(), k);

                for canonical in kmerizer {
                    set.insert(canonical);
                }
            }
        }

        Self { set, k }
    }

    #[cfg(feature = "parallel")]
    pub fn from_fastq<R>(input: R, k: u8) -> Self
    where
        R: std::io::BufRead,
    {
        let mut set = rustc_hash::FxHashSet::default();

        let mut reader = noodles::fastq::Reader::new(input);
        let mut iter = reader.records();
        let mut records = Vec::with_capacity(8192);

        let mut end = true;
        while end {
            log::info!("Start populate buffer");
            end = populate_bufferq(&mut iter, &mut records, 8192);
            log::info!("End populate buffer {}", records.len());

            set.extend(
                records
                    .par_iter()
                    .filter(|record| record.sequence().len() >= k as usize)
                    .map(|record| {
                        let mut set = rustc_hash::FxHashSet::default();
                        let kmerizer =
                            cocktail::tokenizer::Canonical::new(record.sequence().as_ref(), k);
                        for canonical in kmerizer {
                            set.insert(canonical);
                        }
                        set
                    })
                    .reduce(
                        || rustc_hash::FxHashSet::default(),
                        |mut x, y| {
                            x.extend(y.iter());
                            x
                        },
                    ),
            );
        }

        Self { set, k }
    }
}

#[cfg(feature = "parallel")]
/// Populate record buffer with content of iterator
fn populate_buffer<R>(
    iter: &mut noodles::fasta::reader::Records<'_, R>,
    records: &mut Vec<noodles::fasta::Record>,
    record_buffer: u64,
) -> bool
where
    R: std::io::BufRead,
{
    records.clear();

    for i in 0..record_buffer {
        if let Some(Ok(record)) = iter.next() {
            records.push(record);
        } else {
            records.truncate(i as usize);
            return false;
        }
    }

    true
}

#[cfg(feature = "parallel")]
/// Populate record buffer with content of iterator
fn populate_bufferq<R>(
    iter: &mut noodles::fastq::reader::Records<'_, R>,
    records: &mut Vec<noodles::fastq::Record>,
    record_buffer: u64,
) -> bool
where
    R: std::io::BufRead,
{
    records.clear();

    for i in 0..record_buffer {
        if let Some(Ok(record)) = iter.next() {
            records.push(record);
        } else {
            records.truncate(i as usize);
            return false;
        }
    }

    true
}

impl set::KmerSet for Hash {
    fn get(&self, kmer: u64) -> bool {
        self.set.contains(&cocktail::kmer::canonical(kmer, self.k))
    }

    fn k(&self) -> u8 {
        self.k
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    static FILE: &[u8] = b">1\nACGTGGGAATTGTGGCCACATCACGAGGTCCTGCGTATTGACGACTGTAAAGCGAGTGGCCGTGGAATTTCAAGCTCAATTAGCCGAACCAATCCGCCTA";

    #[test]
    fn canonical() {
        let file = std::io::Cursor::new(FILE);

        let hash = Hash::from_fasta(file, 11);

        let set: crate::set::BoxKmerSet = Box::new(hash);

        let mut records = bio::io::fasta::Reader::new(FILE).records();
        for cano in cocktail::tokenizer::Canonical::new(records.next().unwrap().unwrap().seq(), 11)
        {
            assert!(set.get(cano))
        }
    }

    #[test]
    fn forward() {
        let file = std::io::Cursor::new(FILE);

        let hash = Hash::from_fasta(file, 11);

        let set: crate::set::BoxKmerSet = Box::new(hash);

        let mut records = bio::io::fasta::Reader::new(FILE).records();
        for kmer in cocktail::tokenizer::Tokenizer::new(records.next().unwrap().unwrap().seq(), 11)
        {
            assert!(set.get(kmer))
        }
    }

    #[test]
    fn absence() {
        let file = std::io::Cursor::new(FILE);

        let hash = Hash::from_fasta(file, 11);

        let set: crate::set::BoxKmerSet = Box::new(hash);

        assert!(!set.get(0));
    }

    #[test]
    fn k() {
        let file = std::io::Cursor::new(FILE);

        let hash = Hash::from_fasta(file, 11);

        let set: crate::set::BoxKmerSet = Box::new(hash);

        assert_eq!(set.k(), 11);
    }
}
