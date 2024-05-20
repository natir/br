//! Brutal Rewrite

#![warn(missing_docs)]

/* std use */

/* crates use */

use anyhow::Context as _;
use clap::Parser as _;

/* project use */
use br::cli;
use br::error;
use br::set;

fn main() -> error::Result<()> {
    // parse cli
    let params = cli::Command::parse();

    // Setup logger
    stderrlog::new()
        .module(module_path!())
        .quiet(params.quiet())
        .verbosity(params.verbosity())
        .timestamp(params.timestamp())
        .init()
        .context("stderrlog already create a logger")?;

    #[cfg(feature = "parallel")]
    rayon::ThreadPoolBuilder::new()
        .num_threads(params.threads())
        .build_global()?;

    let kmer_set = match params.subcommand {
        cli::SubCommand::Count(ref subparams) => count(subparams),
        cli::SubCommand::Fasta(ref subparams) => fasta(subparams),
        cli::SubCommand::Solid(ref subparams) => solid(subparams),
        cli::SubCommand::LargeKmer(ref subparams) => large_kmer(subparams),
    }?;

    let methods = br::build_methods(
        params.corrections(),
        &kmer_set,
        params.confirm(),
        params.max_search(),
    );

    br::run_correction(
        &mut params.inputs()?,
        &mut params.outputs()?,
        methods,
        params.two_side(),
        params.record_buffer(),
    )?;

    Ok(())
}

fn count(subparams: &br::cli::Count) -> error::Result<Box<dyn br::set::KmerSet>> {
    let counter =
        pcon::counter::Counter::<pcon::CountTypeNoAtomic>::from_stream(subparams.inputs()?)?;

    count2solid(
        counter.raw(),
        subparams.abundance(),
        subparams.abundance_selection(),
        counter.k(),
    )
}

fn fasta(subparams: &br::cli::Fasta) -> error::Result<Box<dyn br::set::KmerSet>> {
    let mut counter = pcon::counter::Counter::<pcon::CountType>::new(subparams.kmer_size());
    counter.count_fasta(subparams.inputs()?, 8192);

    let raw = unsafe {
        std::mem::transmute::<&[pcon::CountType], &[pcon::CountTypeNoAtomic]>(counter.raw())
    };
    count2solid(
        raw,
        subparams.abundance(),
        subparams.abundance_selection(),
        counter.k(),
    )
}

fn count2solid(
    counter: &[pcon::CountTypeNoAtomic],
    abundance: Option<pcon::CountTypeNoAtomic>,
    abundance_selection: Option<cli::AbundanceSelection>,
    kmer_size: pcon::CountTypeNoAtomic,
) -> error::Result<Box<dyn br::set::KmerSet>> {
    let spectrum = pcon::spectrum::Spectrum::from_count(counter);

    let abundance = match (abundance, abundance_selection) {
        (Some(x), _) => Ok(x),
        (_, Some(cli::AbundanceSelection::FirstMinimum)) => spectrum
            .get_threshold(pcon::spectrum::ThresholdMethod::FirstMinimum, 0.0)
            .ok_or(error::Error::ComputeAbundanceThreshold),
        (_, Some(cli::AbundanceSelection::Rarefaction { percent })) => spectrum
            .get_threshold(pcon::spectrum::ThresholdMethod::Rarefaction, percent)
            .ok_or(error::Error::ComputeAbundanceThreshold),
        (_, Some(cli::AbundanceSelection::PercentLeast { percent })) => spectrum
            .get_threshold(pcon::spectrum::ThresholdMethod::PercentAtLeast, percent)
            .ok_or(error::Error::ComputeAbundanceThreshold),
        (_, Some(cli::AbundanceSelection::PercentMost { percent })) => spectrum
            .get_threshold(pcon::spectrum::ThresholdMethod::PercentAtMost, percent)
            .ok_or(error::Error::ComputeAbundanceThreshold),
        (None, None) => Err(error::Error::AbundanceThresholdOrAbundanceMethod),
    }?;

    Ok(Box::new(set::Pcon::new(pcon::solid::Solid::from_count(
        kmer_size, counter, abundance,
    ))))
}

fn solid(subparams: &br::cli::Solid) -> error::Result<Box<dyn br::set::KmerSet>> {
    let set = match subparams.format() {
        cli::SolidInput::Solid => set::Pcon::from_pcon_solid(subparams.input()?)?,
        cli::SolidInput::Csv => set::Pcon::from_csv(
            subparams.input()?,
            subparams
                .kmer_size()
                .ok_or(error::Error::SolidRequireKmerSize)?,
        )?,
        cli::SolidInput::Fasta => set::Pcon::from_fasta(
            subparams.input()?,
            subparams
                .kmer_size()
                .ok_or(error::Error::SolidRequireKmerSize)?,
        ),
        cli::SolidInput::Fastq => set::Pcon::from_fastq(
            subparams.input()?,
            subparams
                .kmer_size()
                .ok_or(error::Error::SolidRequireKmerSize)?,
        ),
        #[cfg(feature = "kff")]
        cli::SolidInput::Kff => todo!(),
    };

    Ok(Box::new(set))
}

fn large_kmer(subparams: &br::cli::LargeKmer) -> error::Result<Box<dyn br::set::KmerSet>> {
    let set = match subparams.format() {
        cli::LargeKmerInput::Csv => set::Hash::from_csv(subparams.input()?, subparams.kmer_size())?,
        cli::LargeKmerInput::Fasta => {
            set::Hash::from_fasta(subparams.input()?, subparams.kmer_size())
        }
        cli::LargeKmerInput::Fastq => {
            set::Hash::from_fastq(subparams.input()?, subparams.kmer_size())
        }
        #[cfg(feature = "kff")]
        cli::LargeKmerInput::Kff => todo!(),
    };

    Ok(Box::new(set))
}
