//! Command Line Interface declaration of project br

/* std use */
use std::io::Read as _;

/* crate use */

/* project use */
use crate::error;

#[derive(std::clone::Clone, std::fmt::Debug, clap::ValueEnum)]
pub enum CorrectionMethod {
    One,
    Two,
    Graph,
    Greedy,
    GapSize,
}

/// Brutal Rewrite, a kmer read corrector
#[derive(clap::Parser, std::fmt::Debug)]
#[clap(
    version = "0.1",
    author = "Pierre Marijon <pierre.marijon@hhu.de>",
    about = "Br: Brutal rewrite a simple long read corrector based on kmer spectrum methode"
)]
pub struct Command {
    /// Path to inputs, default read stdin
    #[clap(short = 'i', long = "inputs")]
    inputs: Option<Vec<std::path::PathBuf>>,

    /// Path to output, default read stdout
    #[clap(short = 'o', long = "outputs")]
    outputs: Option<Vec<std::path::PathBuf>>,

    /// Correct in two side
    #[clap(short = 's', long = "two-side")]
    two_side: bool,

    /// Correction method
    #[clap(short = 'c', long = "corrections")]
    corrections: Option<Vec<CorrectionMethod>>,

    /// Number of kmer required to validate correction
    #[clap(short = 'C', long = "confirm")]
    confirm: Option<u8>,

    /// Number of base we use to try correct error, default value is '7'
    #[clap(short = 'M', long = "max-search")]
    max_search: Option<u8>,

    /// Number of sequence record load in buffer, default 8192
    #[clap(short = 'b', long = "record_buffer")]
    pub record_buffer: Option<u64>,

    /// SubCommand
    #[clap(subcommand)]
    pub subcommand: SubCommand,

    #[cfg(feature = "parallel")]
    /// Number of theard use 0 use all avaible core, default value 0
    #[clap(short = 't', long = "threads")]
    threads: Option<usize>,

    /// Silence all output
    #[clap(short = 'q', long = "quiet")]
    quiet: bool,

    /// Verbose mode (-v, -vv, -vvv, etc)
    #[clap(short = 'v', long = "verbosity", action = clap::ArgAction::Count)]
    verbosity: u8,

    /// Timestamp (sec, ms, ns, none)
    #[clap(short = 'T', long = "timestamp")]
    ts: Option<stderrlog::Timestamp>,
}

impl Command {
    /// Get inputs
    pub fn inputs(&self) -> error::Result<Vec<Box<dyn std::io::BufRead>>> {
        let mut inputs: Vec<Box<dyn std::io::BufRead>> = vec![];

        match &self.inputs {
            None => Ok(vec![Box::new(std::io::stdin().lock())]),
            Some(paths) => {
                for path in paths {
                    inputs.push(open(path)?);
                }

                Ok(inputs)
            }
        }
    }

    /// Get output
    pub fn outputs(&self) -> error::Result<Vec<Box<dyn std::io::Write>>> {
        let mut outputs: Vec<Box<dyn std::io::Write>> = vec![];

        match &self.outputs {
            None => {
                if outputs.is_empty() {
                    outputs.push(Box::new(std::io::BufWriter::new(std::io::stdout())))
                }
            }
            Some(paths) => {
                for path in paths {
                    outputs.push(create(path)?);
                }
            }
        }

        Ok(outputs)
    }

    /// Get two_side
    pub fn two_side(&self) -> bool {
        self.two_side
    }

    /// Get correction method
    pub fn corrections(&self) -> Vec<CorrectionMethod> {
        match &self.corrections {
            Some(ms) => ms.to_vec(),
            None => vec![
                CorrectionMethod::One,
                CorrectionMethod::Two,
                CorrectionMethod::Graph,
                CorrectionMethod::Greedy,
                CorrectionMethod::GapSize,
            ],
        }
    }

    /// Get confirm
    pub fn confirm(&self) -> u8 {
        self.confirm.unwrap_or(5)
    }

    /// Get confirm
    pub fn max_search(&self) -> u8 {
        self.max_search.unwrap_or(7)
    }

    /// Get record buffer
    pub fn record_buffer(&self) -> u64 {
        self.record_buffer.unwrap_or(8192)
    }

    /// Get number of thread
    #[cfg(feature = "parallel")]
    pub fn threads(&self) -> usize {
        self.threads.unwrap_or(0)
    }

    /// Get verbosity level
    pub fn verbosity(&self) -> usize {
        self.verbosity as usize
    }

    /// Get quiet
    pub fn quiet(&self) -> bool {
        self.quiet
    }

    /// Get timestamp granularity
    pub fn timestamp(&self) -> stderrlog::Timestamp {
        self.ts.unwrap_or(stderrlog::Timestamp::Off)
    }
}

/// Enumeration of subcommand
#[derive(clap::Subcommand, std::fmt::Debug)]
pub enum SubCommand {
    /// With Count
    Count(Count),

    /// With Fasta
    Fasta(Fasta),

    /// With Solid
    Solid(Solid),

    /// Large Kmer mode
    LargeKmer(LargeKmer),
}

/// SubCommand Count
#[derive(clap::Args, std::fmt::Debug)]
pub struct Count {
    /// Path to fasta kmer inputs
    #[clap(short = 'i', long = "inputs")]
    inputs: std::path::PathBuf,

    /// Minimal abundance, default value 0
    #[clap(short = 'a', long = "abundance")]
    abundance: Option<pcon::CountTypeNoAtomic>,

    /// Abundance selection method
    #[clap(subcommand)]
    abundance_selection: Option<AbundanceSelection>,
}

impl Count {
    /// Get inputs
    pub fn inputs(&self) -> error::Result<Box<dyn std::io::BufRead>> {
        let mut handle: Box<dyn std::io::Read> = Box::new(std::io::Cursor::new(vec![]));

        for path in &self.inputs {
            let (file, _compression) = niffler::get_reader(Box::new(std::fs::File::open(path)?))?;
            handle = Box::new(handle.chain(file));
        }

        Ok(Box::new(std::io::BufReader::new(handle)))
    }

    /// Get abundance
    pub fn abundance(&self) -> Option<pcon::CountTypeNoAtomic> {
        self.abundance
    }

    /// Get abundance selection method
    pub fn abundance_selection(&self) -> Option<AbundanceSelection> {
        self.abundance_selection
    }
}

/// Enumeration of abundance selection
#[derive(clap::Subcommand, std::fmt::Debug, std::clone::Clone, std::marker::Copy)]
pub enum AbundanceSelection {
    /// First Minimum
    FirstMinimum,

    /// Rarefaction
    Rarefaction { percent: f64 },

    /// PercentMost
    PercentMost { percent: f64 },

    /// PercentLeast
    PercentLeast { percent: f64 },
}

/// SubCommand Fasta
#[derive(clap::Args, std::fmt::Debug)]
pub struct Fasta {
    /// Path to fasta kmer inputs
    #[clap(short = 'i', long = "inputs")]
    inputs: Vec<std::path::PathBuf>,

    /// Size of kmer
    #[clap(short = 'k', long = "kmer-size")]
    kmer_size: u8,

    /// Minimal abundance, default value 0
    #[clap(short = 'a', long = "abundance")]
    abundance: Option<pcon::CountTypeNoAtomic>,

    /// Abundance selection method
    #[clap(subcommand)]
    pub abundance_selection: Option<AbundanceSelection>,
}

impl Fasta {
    /// Get inputs
    pub fn inputs(&self) -> error::Result<Box<dyn std::io::BufRead>> {
        let mut handle: Box<dyn std::io::Read> = Box::new(std::io::Cursor::new(vec![]));

        for path in &self.inputs {
            let (file, _compression) = niffler::get_reader(Box::new(std::fs::File::open(path)?))?;
            handle = Box::new(handle.chain(file));
        }

        Ok(Box::new(std::io::BufReader::new(handle)))
    }

    /// Get size of kmer
    pub fn kmer_size(&self) -> u8 {
        self.kmer_size - (!(self.kmer_size & 0b1) & 0b1)
    }

    /// Get abundance
    pub fn abundance(&self) -> Option<pcon::CountTypeNoAtomic> {
        self.abundance
    }

    /// Get abundance selection method
    pub fn abundance_selection(&self) -> Option<AbundanceSelection> {
        self.abundance_selection
    }
}

#[derive(clap::ValueEnum, std::clone::Clone, std::fmt::Debug)]
pub enum SolidInput {
    Solid,
    Csv,
    Fasta,
    Fastq,
    #[cfg(feature = "kff")]
    Kff,
}

/// SubCommand Solid
#[derive(clap::Args, std::fmt::Debug)]
pub struct Solid {
    /// Path to input
    #[clap(short = 'i', long = "input")]
    input: std::path::PathBuf,

    /// Input type
    #[clap(short = 'f', long = "format")]
    format: SolidInput,

    /// Size of kmer
    #[clap(short = 'k', long = "kmer-size")]
    kmer_size: Option<u8>,
}

impl Solid {
    /// Get input
    pub fn input(&self) -> error::Result<Box<dyn std::io::BufRead + std::marker::Send>> {
        open_send(&self.input)
    }

    /// Get format
    pub fn format(&self) -> SolidInput {
        self.format.clone()
    }

    /// Get size of kmer
    pub fn kmer_size(&self) -> Option<u8> {
        self.kmer_size
    }
}

#[derive(clap::ValueEnum, std::clone::Clone, std::fmt::Debug)]
pub enum LargeKmerInput {
    Csv,
    Fasta,
    Fastq,
    #[cfg(feature = "kff")]
    Kff,
}

/// SubCommand LargeKmer
#[derive(clap::Args, std::fmt::Debug)]
pub struct LargeKmer {
    /// Path to input
    #[clap(short = 'i', long = "input")]
    input: std::path::PathBuf,

    /// Input type
    #[clap(short = 'f', long = "format")]
    format: LargeKmerInput,

    /// Size of kmer
    #[clap(short = 'k', long = "kmer-size")]
    kmer_size: u8,
}

impl LargeKmer {
    /// Get input
    pub fn input(&self) -> error::Result<Box<dyn std::io::BufRead + std::marker::Send>> {
        match self.format {
            LargeKmerInput::Csv => open_send(&self.input),
            LargeKmerInput::Fasta => open_send(&self.input),
            LargeKmerInput::Fastq => open_send(&self.input),
            #[cfg(feature = "kff")]
            LargeKmerInput::Kff => open_send(&self.input),
        }
    }

    /// Get format
    pub fn format(&self) -> LargeKmerInput {
        self.format.clone()
    }

    /// Get size of kmer
    pub fn kmer_size(&self) -> u8 {
        self.kmer_size
    }
}

fn create<P>(path: P) -> error::Result<Box<dyn std::io::Write + std::marker::Send>>
where
    P: std::convert::AsRef<std::path::Path>,
{
    let file = std::fs::File::create(path)?;
    let buffer = std::io::BufWriter::new(file);
    let boxed = Box::new(buffer);

    Ok(boxed)
}

fn open<P>(path: P) -> error::Result<Box<dyn std::io::BufRead>>
where
    P: std::convert::AsRef<std::path::Path>,
{
    let file = niffler::get_reader(Box::new(std::fs::File::open(path)?))?.0;
    let buffer = std::io::BufReader::new(file);
    let boxed = Box::new(buffer);

    Ok(boxed)
}

fn open_send<P>(path: P) -> error::Result<Box<dyn std::io::BufRead + std::marker::Send>>
where
    P: std::convert::AsRef<std::path::Path>,
{
    let file = niffler::send::get_reader(Box::new(std::fs::File::open(path)?))?.0;
    let buffer = std::io::BufReader::new(file);
    let boxed = Box::new(buffer);

    Ok(boxed)
}

#[cfg(test)]
mod tests {
    /*std use */

    /* project use */
    use super::*;

    #[cfg(not(feature = "parallel"))]
    #[test]
    fn basic() {
        let subcmd = Fasta {
            inputs: vec![],
            kmer_size: 14,
            abundance: Some(2),
            abundance_selection: None,
        };

        let cmd = Command {
            inputs: None,
            outputs: None,
            two_side: true,
            confirm: Some(5),
            max_search: Some(7),
            record_buffer: Some(8192),
            corrections: None,
            verbosity: 3,
            quiet: false,
            ts: None,
            subcommand: SubCommand::Fasta(subcmd),
        };

        assert_eq!(cmd.verbosity(), 3);
        assert!(!cmd.quiet());
        assert!(matches!(cmd.timestamp(), stderrlog::Timestamp::Off));

        match cmd.subcommand {
            SubCommand::Fasta(subcmd) => {
                assert_eq!(subcmd.kmer_size(), 13);
                assert_eq!(subcmd.abundance(), Some(2));
                assert!(matches!(subcmd.abundance_selection(), None));
            }
            _ => unreachable!(),
        }
    }

    #[cfg(feature = "parallel")]
    #[test]
    fn basic_parallel() {
        let subcmd = Fasta {
            inputs: vec![],
            kmer_size: 14,
            abundance: Some(2),
            abundance_selection: None,
        };

        let cmd = Command {
            inputs: None,
            outputs: None,
            two_side: true,
            confirm: Some(5),
            max_search: Some(7),
            record_buffer: Some(8192),
            corrections: None,
            verbosity: 3,
            quiet: false,
            ts: None,
            threads: Some(8),
            subcommand: SubCommand::Fasta(subcmd),
        };

        assert_eq!(cmd.verbosity(), 3);
        assert!(!cmd.quiet());
        assert!(matches!(cmd.timestamp(), stderrlog::Timestamp::Off));
        assert_eq!(cmd.threads(), 8);

        match cmd.subcommand {
            SubCommand::Fasta(subcmd) => {
                assert_eq!(subcmd.kmer_size(), 13);
                assert_eq!(subcmd.abundance(), Some(2));
                assert!(matches!(subcmd.abundance_selection(), None));
            }
            _ => unreachable!(),
        }
    }
}
