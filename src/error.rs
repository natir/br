/*
Copyright (c) 2020 Pierre Marijon <pmarijon@mpi-inf.mpg.de>

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
use thiserror::Error;

/// All error produce by Pcon
#[derive(Debug, Error)]
pub enum Error {
    /// See enum [Cli]
    #[error(transparent)]
    Cli(#[from] Cli),

    /// See enum [IO]
    #[error(transparent)]
    IO(#[from] IO),

    #[error("Can't compute minimal abundance")]
    CantComputeAbundance,
}

/// Error emmit durring Cli parsing
#[derive(Debug, Error)]
pub enum Cli {
    /// Number of inputs and outputs must be the same
    #[error("Kmer size must be odd")]
    NotSameNumberOfInAndOut,

    #[error("You must provide a solidity path '-s', a kmer solid path '-S' or a kmer length '-k'")]
    NoSolidityNoKmer,

    #[error("If you provide kmer solid path '-S' you must provide a kmer length '-k'")]
    KmerSolidNeedK,

    #[error("Abundance method threshold can't be parse")]
    CantParseAbundanceMethod,
}

/// Error emmit when pcon try to work with file
#[repr(C)]
#[derive(Debug, Error)]
pub enum IO {
    /// We can't create file. In C binding it's equal to 0
    #[error("We can't create file")]
    CantCreateFile,

    /// We can't open file. In C binding it's equal to 1
    #[error("We can't open file")]
    CantOpenFile,

    /// Error durring write in file. In C binding it's equal to 2
    #[error("Error durring write")]
    ErrorDurringWrite,

    /// Error durring read file. In C binding it's equal to 3
    #[error("Error durring read")]
    ErrorDurringRead,

    /// No error, this exist only for C binding it's the value of a new error pointer
    #[error("Isn't error if you see this please contact the author with this message and a description of what you do with pcon")]
    NoError,
}
