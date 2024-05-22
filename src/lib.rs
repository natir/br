//! Brutal rewrite

/* std use */

/* crates use */

/* project use */

/* mod declaration */
pub mod cli;
pub mod correct;
pub mod error;
pub mod set;

/* crate use */
#[cfg(feature = "parallel")]
use rayon::prelude::*;

/* local use */

#[cfg(not(feature = "parallel"))]
pub fn run_correction<'a>(
    inputs: &mut [Box<dyn std::io::BufRead>],
    outputs: &mut [Box<dyn std::io::Write>],
    methods: Vec<Box<dyn correct::Corrector + Sync + Send + 'a>>,
    two_side: bool,
    _record_buffer_len: u64,
) -> error::Result<()> {
    for (input, output) in inputs.iter_mut().zip(outputs.iter_mut()) {
        let mut reader = noodles::fasta::Reader::new(input);
        let mut writer = noodles::fasta::Writer::new(output);

        let mut records = reader.records();

        while let Some(Ok(record)) = records.next() {
            log::debug!(
                "begin correct read {} {}",
                String::from_utf8(record.name().to_vec()).unwrap(),
                record.sequence().len()
            );
            let seq = record.sequence();

            let mut correct = seq.as_ref().to_vec();
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

            writer.write_record(&noodles::fasta::Record::new(
                record.definition().clone(),
                correct.into(),
            ))?;
            log::debug!(
                "end correct read {}",
                String::from_utf8(record.name().to_vec()).unwrap()
            );
        }
    }

    Ok(())
}

#[cfg(feature = "parallel")]
pub fn run_correction<'a>(
    inputs: &mut [Box<dyn std::io::BufRead>],
    outputs: &mut [Box<dyn std::io::Write>],
    methods: Vec<Box<dyn correct::Corrector + Sync + Send + 'a>>,
    two_side: bool,
    record_buffer_len: u64,
) -> error::Result<()> {
    for (input, output) in inputs.iter_mut().zip(outputs.iter_mut()) {
        let mut reader = noodles::fasta::Reader::new(input);
        let mut writer = noodles::fasta::Writer::new(output);

        let mut iter = reader.records();
        let mut records = Vec::with_capacity(record_buffer_len as usize);
        let mut corrected: Vec<error::Result<noodles::fasta::Record>>;

        let mut end = false;
        while !end {
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
                    log::debug!(
                        "begin correct read {} {}",
                        String::from_utf8(record.name().to_vec()).unwrap(),
                        record.sequence().len()
                    );

                    let seq = record.sequence();

                    let mut correct = seq.as_ref().to_vec();
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

                    log::debug!(
                        "end correct read {}",
                        String::from_utf8(record.name().to_vec()).unwrap()
                    );
                    Ok(noodles::fasta::Record::new(
                        record.definition().clone(),
                        correct.into(),
                    ))
                })
                .collect();

            for corr in corrected {
                writer.write_record(&corr?)?
            }

            records.clear();
        }
    }

    Ok(())
}

pub fn build_methods<'a>(
    params: Vec<cli::CorrectionMethod>,
    solid: &'a set::BoxKmerSet,
    confirm: u8,
    max_search: u8,
) -> Vec<Box<dyn correct::Corrector + Sync + Send + 'a>> {
    let mut methods: Vec<Box<dyn correct::Corrector + Sync + Send + 'a>> = Vec::new();

    for method in params {
        match method {
            cli::CorrectionMethod::One => methods.push(Box::new(correct::One::new(solid, confirm))),
            cli::CorrectionMethod::Two => methods.push(Box::new(correct::Two::new(solid, confirm))),
            cli::CorrectionMethod::Graph => methods.push(Box::new(correct::Graph::new(solid))),
            cli::CorrectionMethod::Greedy => {
                methods.push(Box::new(correct::Greedy::new(solid, max_search, confirm)))
            }
            cli::CorrectionMethod::GapSize => {
                methods.push(Box::new(correct::GapSize::new(solid, confirm)))
            }
        }
    }

    methods
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn methods_list() {
        // Not perfect test

        let set: set::BoxKmerSet = Box::new(set::Pcon::new(pcon::solid::Solid::new(5)));

        let methods = build_methods(
            vec![cli::CorrectionMethod::One, cli::CorrectionMethod::Two],
            &set,
            2,
            5,
        );

        assert_eq!(methods.len(), 2);

        let methods = build_methods(
            vec![
                cli::CorrectionMethod::One,
                cli::CorrectionMethod::Two,
                cli::CorrectionMethod::Graph,
                cli::CorrectionMethod::Greedy,
                cli::CorrectionMethod::GapSize,
                cli::CorrectionMethod::GapSize,
            ],
            &set,
            2,
            5,
        );

        assert_eq!(methods.len(), 6);
    }
}
