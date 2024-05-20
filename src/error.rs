//! Error struct of project br

/* std use */

/* crate use */
use anyhow;
use thiserror;

/* project use */

/// Enum to manage error
#[derive(std::fmt::Debug, thiserror::Error)]
pub enum Error {
    /// Error in logging system configuration
    #[error(transparent)]
    Log(#[from] log::SetLoggerError),

    /// Error in rayon thread pool build
    #[cfg(feature = "parallel")]
    #[error(transparent)]
    RayonThreadPool(#[from] rayon::ThreadPoolBuildError),

    /// Error in input output
    #[error(transparent)]
    IO(#[from] std::io::Error),

    /// Csv didn't contains first column
    #[error("Csv input not contains first column")]
    CsvMissingFirstColumn,

    /// Found minimal threshold failled
    #[error("Br can't compute abundance threshold choose another method")]
    ComputeAbundanceThreshold,

    /// In count and reads subcommand user should set minimum abundance or abundance selection method
    #[error("In count and reads subcommand user should set minimum abundance or abundance selection method")]
    AbundanceThresholdOrAbundanceMethod,

    /// In solid mode csv, fasta and fastq format require kmer size
    #[error("In solid mode csv, fasta and fastq format require kmer size")]
    SolidRequireKmerSize,
}

/// Alias of result
pub type Result<T> = anyhow::Result<T>;
