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
use anyhow::{anyhow, Result};
use clap::Clap;

use br::error::Cli::*;
use br::error::*;
use br::*;

#[derive(Clap, Debug)]
#[clap(
    version = "0.1",
    author = "Pierre Marijon <pierre.marijon@hhu.de>",
    about = "Br: Brutal rewrite a simple long read corrector based on kmer spectrum methode"
)]
pub struct Command {
    #[clap(
        short = 's',
        long = "solidity",
        about = "solidity bitfield produce by pcon"
    )]
    pub solidity: Option<String>,

    #[clap(
        short = 'k',
        long = "kmer",
        about = "kmer length if you didn't provide solidity path you must give a kmer length"
    )]
    pub kmer_size: Option<u8>,

    #[clap(short = 'i', long = "inputs", about = "fasta file to be correct")]
    pub inputs: Vec<String>,

    #[clap(
        short = 'o',
        long = "outputs",
        about = "path where corrected read was write"
    )]
    pub outputs: Vec<String>,

    #[clap(
        short = 'a',
        long = "abundance",
        about = "if you want choose the minimum abundance you can set this parameter"
    )]
    pub abundance: Option<u8>,

    #[clap(
	short = 'm',
	long = "method",
	possible_values = &["one", "graph", "greedy", "gap_size"],
	about = "correction method used, methods are applied in the order you specify, default value is 'one'"
    )]
    pub methods: Option<Vec<String>>,

    #[clap(
        short = 'c',
        long = "confirm",
        about = "number of kmer need to be solid after one, greedy correction to validate it, default value is '2'"
    )]
    pub confirm: Option<u8>,

    #[clap(
        short = 'M',
        long = "max-search",
        about = "number of base we use to try correct error, default value is '7'"
    )]
    pub max_search: Option<u8>,

    #[clap(
        short = 'n',
        long = "not-two-side",
        about = "if this flag is set br correct only in forward orientation"
    )]
    pub two_side: bool,

    #[clap(
        short = 't',
        long = "threads",
        about = "Number of thread use by br, 0 use all avaible core, default value 0"
    )]
    pub threads: Option<usize>,

    #[clap(
        short = 'b',
        long = "record_buffer",
        about = "Number of sequence record load in buffer, default 8192"
    )]
    pub record_buffer: Option<usize>,

    #[clap(
        short = 'v',
        long = "verbosity",
        parse(from_occurrences),
        about = "verbosity level also control by environment variable BR_LOG if flag is set BR_LOG value is ignored"
    )]
    pub verbosity: i8,
}

fn main() -> Result<()> {
    let params = Command::parse();

    if params.inputs.len() != params.outputs.len() {
        return Err(anyhow!(Error::Cli(NotSameNumberOfInAndOut)));
    }

    if let Some(level) = cli::i82level(params.verbosity) {
        env_logger::builder()
            .format_timestamp(Some(env_logger::fmt::TimestampPrecision::Millis))
            .filter_level(level.to_level_filter())
            .init();
    } else {
        env_logger::Builder::from_env("BR_LOG")
            .format_timestamp(Some(env_logger::fmt::TimestampPrecision::Millis))
            .init();
    }

    let confirm = params.confirm.unwrap_or(2);
    let max_search = params.max_search.unwrap_or(7);
    let record_buffer = params.record_buffer.unwrap_or(8192);

    if let Some(threads) = params.threads {
        log::info!("Set number of threads to {}", threads);

        set_nb_threads(threads);
    }

    let solid = cli::read_or_compute_solidity(
        params.solidity,
        params.kmer_size,
        &params.inputs,
        record_buffer,
        params.abundance,
    )?;

    let methods = br::build_methods(params.methods, &solid, confirm, max_search);

    br::run_correction(
        &params.inputs,
        &params.outputs,
        methods,
        params.two_side,
        record_buffer,
    )?;

    Ok(())
}
