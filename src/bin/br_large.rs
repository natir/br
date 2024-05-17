/*
Copyright (c) 2021 Pierre Marijon <pierre.marijon@hhu.de>

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
use clap::Parser;

use br::error::IO::*;
use br::error::*;
use br::*;

#[derive(clap::Parser, Debug)]
#[clap(
    version = "0.1",
    author = "Pierre Marijon <pierre.marijon@hhu.de>",
    about = "A version of br where memory usage depend on number of kmer not kmer length"
)]
pub struct Command {
    /// fasta file to be correct
    #[clap(short = 'i', long = "inputs")]
    pub inputs: Vec<String>,

    /// path where corrected read was write
    #[clap(short = 'o', long = "outputs")]
    pub outputs: Vec<String>,

    /// use kmer present in fasta file as solid kmer and store them in HashSet
    #[clap(short = 'S', long = "kmer-solid")]
    pub kmer_solid: Vec<String>,

    /// kmer length lower or equal to 32
    #[clap(short = 'k', long = "kmer")]
    pub kmer_size: Option<u8>,

    /// correction method used, methods are applied in the order you specify, default value is 'one'
    #[clap(short = 'm', long = "method")]
    pub methods: Option<Vec<br::CorrectionMethod>>,

    /// number of kmer need to be solid after one, greedy correction to validate it, default value is '2'
    #[clap(short = 'c', long = "confirm")]
    pub confirm: Option<u8>,

    /// number of base we use to try correct error, default value is '7'
    #[clap(short = 'M', long = "max-search")]
    pub max_search: Option<u8>,

    /// if this flag is set br correct only in forward orientation
    #[clap(short = 'n', long = "not-two-side")]
    pub two_side: bool,

    /// Number of thread use by br, 0 use all avaible core, default value 0
    #[clap(short = 't', long = "threads")]
    pub threads: Option<usize>,

    /// Number of sequence record load in buffer, default 8192
    #[clap(short = 'b', long = "record_buffer")]
    pub record_buffer: Option<u64>,

    /// verbosity level also control by environment variable BR_LOG if flag is set BR_LOG value is ignored
    #[clap(short = 'v', long = "verbosity", action =  clap::ArgAction::Count)]
    pub verbosity: u8,
}

#[cfg(not(tarpaulin_include))]
fn main() -> Result<()> {
    let params = Command::parse();

    let mut files = Vec::new();

    for path in params.kmer_solid {
        files.push(std::io::BufReader::new(
            std::fs::File::open(&path)
                .with_context(|| Error::IO(CantOpenFile))
                .with_context(|| anyhow!("File {:?}", path.clone()))?,
        ));
    }

    let solid: set::BoxKmerSet = Box::new(set::Hash::new(
        files,
        params.kmer_size.ok_or(Error::Cli(Cli::KmerSolidNeedK))?,
    ));

    let confirm = params.confirm.unwrap_or(2);
    let max_search = params.max_search.unwrap_or(7);
    let record_buffer = params.record_buffer.unwrap_or(8192);

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
