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
use structopt::StructOpt;

#[derive(StructOpt, Debug)]
#[structopt(
    version = "0.1",
    author = "Pierre Marijon <pmarijon@mpi-inf.mpg.de>",
    about = "Br: Brutal rewrite a simple long read corrector based on kmer spectrum methode"
)]
pub struct Command {
    #[structopt(short = "s", long = "solidity", help = "solidity bitfield produce by pcon")]
    pub solidity: String,

    #[structopt(short = "i", long = "inputs", help = "fasta file to be correct")]
    pub inputs: Vec<String>,

    #[structopt(short = "o", long = "outputs", help = "path where corrected read was write")]
    pub outputs: Vec<String>,

    #[structopt(short = "c", long = "confirm", help = "number of kmer need to be solid after correction to validate it")]
    pub confirm: u8,
}

use crate::error::{Cli, Error};
use Cli::*;

pub fn check_params(params: Command) -> Result<Command, Error> {
    if params.inputs.len() != params.outputs.len() {
	return Err(Error::Cli(NotSameNumberOfInAndOut));
    }

    Ok(params)
}
