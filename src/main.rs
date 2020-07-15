/*
Copyright (c) 2019 Pierre Marijon <pmarijon@mmci.uni-saarland.de>

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
use clap::Clap;

use log::info;

use br::error::IO::*;
use br::error::*;
use br::*;

fn main() -> Result<()> {
    let mut params = cli::Command::parse();

    params = cli::check_params(params)?;

    if let Some(level) = cli::i82level(params.verbosity) {
        env_logger::builder()
            .format_timestamp(None)
            .filter_level(level.to_level_filter())
            .init();
    } else {
        env_logger::Builder::from_env("BR_LOG")
            .format_timestamp(None)
            .init();
    }

    let solidity_reader = std::io::BufReader::new(
        std::fs::File::open(&params.solidity)
            .with_context(|| Error::IO(CantOpenFile))
            .with_context(|| anyhow!("File {}", params.solidity.clone()))?,
    );
    let solid = pcon::solid::Solid::deserialize(solidity_reader)?;

    let confirm = if let Some(val) = params.confirm {
        val
    } else {
        2
    };

    for (input, output) in params.inputs.iter().zip(params.outputs) {
        info!("Read file {} write in {}", input, output);

        let reader = bio::io::fasta::Reader::new(std::io::BufReader::new(
            std::fs::File::open(input)
                .with_context(|| Error::IO(CantOpenFile))
                .with_context(|| anyhow!("File {}", input.clone()))?,
        ));

        let mut write = bio::io::fasta::Writer::new(std::io::BufWriter::new(
            std::fs::File::create(&output)
                .with_context(|| Error::IO(CantCreateFile))
                .with_context(|| anyhow!("File {}", output.clone()))?,
        ));

        let mut records = reader.records();
        while let Some(Ok(record)) = records.next() {
            info!("correct read {}", record.id());

            let seq = record.seq();

            let correct = if let Some(methods) = &params.methods {
                let mut tmp = seq.to_vec();

                if methods.contains(&"greedy".to_string()) {
                    tmp = correct::greedy::correct(tmp.as_slice(), &solid, confirm);
                }

                if methods.contains(&"graph".to_string()) {
                    tmp = correct::graph::correct(tmp.as_slice(), &solid);
                }

                if methods.contains(&"gap_size".to_string()) {
                    tmp = correct::gap_size::correct(tmp.as_slice(), &solid, confirm);
                }

                tmp
            } else {
                correct::gap_size::correct(seq, &solid, confirm)
            };

            write
                .write_record(&bio::io::fasta::Record::with_attrs(
                    record.id(),
                    record.desc(),
                    &correct,
                ))
                .with_context(|| Error::IO(IO::ErrorDurringWrite))
                .with_context(|| anyhow!("File {}", input.clone()))?;
        }
    }

    Ok(())
}
