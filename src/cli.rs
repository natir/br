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
use anyhow::{anyhow, Context, Result};
use log::Level;

use crate::error::*;
use crate::set;
use Cli::*;
use IO::*;

pub fn i82level(level: i8) -> Option<Level> {
    match level {
        std::i8::MIN..=0 => None,
        1 => Some(log::Level::Error),
        2 => Some(log::Level::Warn),
        3 => Some(log::Level::Info),
        4 => Some(log::Level::Debug),
        5..=std::i8::MAX => Some(log::Level::Trace),
    }
}

pub fn read_or_compute_solidity(
    solidity_path: Option<String>,
    kmer_size: Option<u8>,
    inputs: &[String],
    record_buffer_len: usize,
    abundance: Option<u8>,
) -> Result<set::BoxKmerSet> {
    if let Some(solidity_path) = solidity_path {
        let solidity_reader = std::io::BufReader::new(
            std::fs::File::open(&solidity_path)
                .with_context(|| Error::IO(CantOpenFile))
                .with_context(|| anyhow!("File {:?}", solidity_path.clone()))?,
        );

        log::info!("Load solidity file");
        Ok(Box::new(set::Pcon::new(pcon::solid::Solid::deserialize(
            solidity_reader,
        )?)))
    } else if let Some(kmer_size) = kmer_size {
        let mut counter = pcon::counter::Counter::new(kmer_size);

        log::info!("Start count kmer from input");
        for input in inputs {
            let fasta = std::io::BufReader::new(
                std::fs::File::open(&input)
                    .with_context(|| Error::IO(CantOpenFile))
                    .with_context(|| anyhow!("File {:?}", input.clone()))?,
            );

            counter.count_fasta(fasta, record_buffer_len);
        }
        log::info!("End count kmer from input");

        let spectrum = pcon::spectrum::Spectrum::from_counter(&counter);

        let abun = if let Some(a) = abundance {
            a
        } else {
            let abundance = spectrum
                .get_threshold(pcon::spectrum::ThresholdMethod::FirstMinimum, 0.0)
                .ok_or(Error::CantComputeAbundance)?;

            spectrum
                .write_histogram(std::io::stdout(), Some(abundance))
                .with_context(|| anyhow!("Error durring write of kmer histograme"))?;
            println!("If this curve seems bad or minimum abundance choose (marked by *) not apopriate set parameter -a");

            abundance as u8
        };

        Ok(Box::new(set::Pcon::new(pcon::solid::Solid::from_counter(
            &counter, abun,
        ))))
    } else {
        Err(Error::Cli(NoSolidityNoKmer).into())
    }
}
