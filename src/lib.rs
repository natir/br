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

/* local mod */
pub mod cli;
pub mod correct;
pub mod error;
pub mod set;

/* crate use */
use anyhow::{anyhow, Context, Result};
use rayon::iter::ParallelBridge;
use rayon::prelude::*;

/* local use */
use error::IO::*;
use error::*;

pub fn run_correction<'a>(
    inputs: &[String],
    outputs: &[String],
    methods: Vec<Box<dyn correct::Corrector + Sync + Send + 'a>>,
    two_side: bool,
    record_buffer_len: usize,
) -> Result<()> {
    for (input, output) in inputs.iter().zip(outputs) {
        log::info!("Read file {} write in {}", input, output);

        let reader = bio::io::fasta::Reader::new(std::io::BufReader::new(
            std::fs::File::open(input)
                .with_context(|| error::Error::IO(CantOpenFile))
                .with_context(|| anyhow!("File {}", input.clone()))?,
        ));

        let mut write = bio::io::fasta::Writer::new(std::io::BufWriter::new(
            std::fs::File::create(&output)
                .with_context(|| error::Error::IO(CantCreateFile))
                .with_context(|| anyhow!("File {}", output.clone()))?,
        ));

        let mut iter = reader.records();
        let mut records = Vec::with_capacity(record_buffer_len);
        let mut corrected: Vec<bio::io::fasta::Record>;

        let mut end = false;
        loop {
            for _ in 0..record_buffer_len {
                if let Some(Ok(record)) = iter.next() {
                    records.push(record);
                } else {
                    end = true;
                    break;
                }
            }

            log::info!("Buffer len: {}", records.len());

            corrected = records
                .drain(..)
                .par_bridge()
                .map(|record| {
                    log::debug!("begin correct read {} {}", record.id(), record.seq().len());

                    let seq = record.seq();

                    let mut correct = seq.to_vec();
                    methods
                        .iter()
                        .for_each(|x| correct = x.correct(correct.as_slice()));

                    if !two_side {
                        correct.reverse();
                        methods
                            .iter()
                            .for_each(|x| correct = x.correct(correct.as_slice()));

                        correct.reverse();
                    }

                    log::debug!("end correct read {}", record.id());
                    bio::io::fasta::Record::with_attrs(record.id(), record.desc(), &correct)
                })
                .collect();

            for corr in corrected {
                write
                    .write_record(&corr)
                    .with_context(|| Error::IO(IO::ErrorDurringWrite))
                    .with_context(|| anyhow!("File {}", output.clone()))?
            }

            records.clear();

            if end {
                break;
            }
        }
    }

    Ok(())
}

pub fn build_methods<'a>(
    params: Option<Vec<String>>,
    solid: &'a set::BoxKmerSet,
    confirm: u8,
    max_search: u8,
) -> Vec<Box<dyn correct::Corrector + Sync + Send + 'a>> {
    let mut methods: Vec<Box<dyn correct::Corrector + Sync + Send + 'a>> = Vec::new();

    if let Some(ms) = params {
        for method in ms {
            match &method[..] {
                "one" => methods.push(Box::new(correct::One::new(solid, confirm))),
                "two" => methods.push(Box::new(correct::Two::new(solid, confirm))),
                "graph" => methods.push(Box::new(correct::Graph::new(&solid))),
                "greedy" => {
                    methods.push(Box::new(correct::Greedy::new(&solid, max_search, confirm)))
                }
                "gap_size" => methods.push(Box::new(correct::GapSize::new(&solid, confirm))),
                _ => unreachable!(),
            }
        }
    } else {
        methods.push(Box::new(correct::One::new(&solid, confirm)));
    }

    methods
}

/// Set the number of threads use by count step
pub fn set_nb_threads(nb_threads: usize) {
    rayon::ThreadPoolBuilder::new()
        .num_threads(nb_threads)
        .build_global()
        .unwrap();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn change_number_of_thread() {
        set_nb_threads(16);
        assert_eq!(rayon::current_num_threads(), 16);
    }
}
