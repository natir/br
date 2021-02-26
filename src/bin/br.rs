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
	possible_values = &["one", "two", "graph", "greedy", "gap_size"],
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

#[cfg(not(tarpaulin_include))]
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

    let solid = if let Some(path) = params.solidity {
        read_solidity(path)?
    } else if let Some(kmer_size) = params.kmer_size {
        let count = count_kmer(&params.inputs, kmer_size, record_buffer)?;

        let threshold = params
            .abundance
            .unwrap_or(compute_abundance_threshold(&count, std::io::stdout())?);

        Box::new(set::Pcon::new(pcon::solid::Solid::from_counter(
            &count, threshold,
        )))
    } else {
        anyhow::bail!(Error::Cli(NoSolidityNoKmer));
    };

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

fn read_solidity<'a>(path: String) -> Result<set::BoxKmerSet<'a>> {
    let solidity_reader = std::io::BufReader::new(
        std::fs::File::open(&path)
            .with_context(|| Error::IO(IO::CantOpenFile))
            .with_context(|| anyhow!("File {:?}", path.clone()))?,
    );

    log::info!("Load solidity file");

    Ok(Box::new(set::Pcon::new(pcon::solid::Solid::deserialize(
        solidity_reader,
    )?)))
}

fn count_kmer(
    inputs: &[String],
    kmer_size: u8,
    record_buffer_len: usize,
) -> Result<pcon::counter::Counter> {
    let mut counter = pcon::counter::Counter::new(kmer_size);

    log::info!("Start count kmer from input");
    for input in inputs {
        let fasta = std::io::BufReader::new(
            std::fs::File::open(&input)
                .with_context(|| Error::IO(IO::CantOpenFile))
                .with_context(|| anyhow!("File {:?}", input.clone()))?,
        );

        counter.count_fasta(fasta, record_buffer_len);
    }
    log::info!("End count kmer from input");

    Ok(counter)
}

fn compute_abundance_threshold<W>(count: &pcon::counter::Counter, out: W) -> Result<u8>
where
    W: std::io::Write,
{
    let spectrum = pcon::spectrum::Spectrum::from_counter(&count);

    let abundance = spectrum
        .get_threshold(pcon::spectrum::ThresholdMethod::FirstMinimum, 0.0)
        .ok_or(Error::CantComputeAbundance)?;

    spectrum
        .write_histogram(out, Some(abundance))
        .with_context(|| anyhow!("Error durring write of kmer histograme"))?;
    println!("If this curve seems bad or minimum abundance choose (marked by *) not apopriate set parameter -a");

    Ok(abundance as u8)
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::io::Seek;
    use std::io::Write;

    use rand::seq::SliceRandom;
    use rand::{Rng, SeedableRng};

    fn init() {
        let _ = env_logger::builder()
            .is_test(true)
            .filter_level(log::LevelFilter::Trace)
            .try_init();
    }

    fn generate_seq_file() -> tempfile::NamedTempFile {
        let mut rng = rand::rngs::SmallRng::seed_from_u64(42);

        let nucs = [b'A', b'C', b'T', b'G'];

        let mut sequence = (0..1000)
            .map(|_| *nucs.choose(&mut rng).unwrap())
            .collect::<Vec<u8>>();

        let mut tmpfile = tempfile::NamedTempFile::new().unwrap();

        writeln!(
            tmpfile.as_file_mut(),
            ">42\n{}",
            std::str::from_utf8(&sequence).unwrap()
        )
        .unwrap();

        for _ in 0..25 {
            for nuc in sequence.iter_mut() {
                if rng.gen_ratio(2, 100) {
                    *nuc = *nucs.choose(&mut rng).unwrap();
                }
            }

            writeln!(
                tmpfile.as_file_mut(),
                ">42\n{}",
                std::str::from_utf8(&sequence).unwrap()
            )
            .unwrap();
        }

        tmpfile.seek(std::io::SeekFrom::Start(0)).unwrap();

        tmpfile
    }

    static COUNT: &[u8] = &[
        5, 31, 139, 8, 0, 0, 0, 0, 0, 4, 255, 5, 192, 88, 111, 210, 96, 140, 18, 100, 144, 140, 66,
        57, 202, 213, 210, 227, 251, 122, 209, 155, 149, 210, 99, 20, 10, 45, 45, 163, 21, 156,
        107, 68, 221, 3, 137, 17, 141, 154, 232, 18, 227, 139, 15, 186, 159, 224, 15, 94, 78, 76,
        83, 27, 34, 168, 118, 148, 79, 51, 105, 75, 153, 164, 57, 61, 52, 44, 92, 198, 11, 32, 67,
        54, 96, 145, 102, 241, 117, 12, 49, 47, 27, 207, 213, 54, 64, 2, 180, 143, 57, 9, 161, 101,
        12, 185, 36, 251, 126, 44, 132, 165, 125, 114, 154, 23, 35, 121, 156, 139, 119, 204, 206,
        161, 217, 58, 42, 106, 113, 224, 81, 128, 27, 6, 198, 78, 233, 104, 93, 85, 187, 130, 198,
        101, 7, 83, 206, 191, 7, 130, 237, 249, 135, 248, 43, 106, 209, 84, 52, 178, 239, 37, 164,
        94, 133, 67, 14, 110, 244, 91, 215, 151, 194, 158, 142, 101, 31, 121, 205, 4, 177, 95, 123,
        38, 187, 235, 233, 213, 146, 41, 81, 119, 52, 93, 236, 202, 35, 95, 82, 192, 246, 115, 201,
        129, 108, 211, 34, 114, 5, 95, 4, 95, 8, 134, 104, 93, 255, 85, 10, 196, 46, 179, 204, 159,
        124, 218, 139, 46, 203, 167, 84, 177, 81, 184, 34, 82, 187, 133, 114, 250, 58, 215, 185,
        176, 225, 30, 132, 237, 74, 125, 61, 17, 106, 78, 94, 55, 230, 76, 210, 234, 244, 90, 239,
        101, 184, 226, 169, 202, 122, 207, 151, 47, 194, 28, 74, 181, 217, 219, 51, 84, 94, 148,
        71, 84, 51, 237, 138, 228, 66, 117, 143, 193, 12, 101, 195, 31, 27, 238, 67, 210, 150, 244,
        181, 84, 0, 145, 165, 0, 106, 248, 28, 151, 31, 67, 27, 161, 176, 229, 125, 247, 92, 5, 47,
        33, 10, 63, 221, 20, 4, 70, 129, 42, 24, 213, 67, 202, 104, 176, 130, 216, 110, 186, 192,
        36, 197, 71, 250, 155, 77, 31, 136, 65, 95, 216, 170, 184, 96, 137, 147, 30, 149, 64, 197,
        5, 209, 12, 5, 222, 120, 96, 48, 209, 196, 251, 142, 201, 176, 103, 139, 165, 7, 16, 236,
        73, 133, 23, 206, 189, 170, 180, 136, 56, 220, 159, 80, 139, 138, 136, 195, 81, 251, 159,
        136, 179, 197, 119, 181, 109, 243, 6, 13, 200, 56, 77, 228, 137, 240, 43, 57, 42, 164, 19,
        14, 5, 246, 149, 250, 230, 182, 242, 159, 43, 136, 23, 148, 125, 242, 178, 181, 96, 84, 56,
        35, 26, 253, 241, 219, 77, 19, 54, 216, 250, 74, 127, 8, 92, 90, 242, 49, 196, 222, 243,
        67, 159, 11, 162, 220, 157, 134, 53, 220, 65, 216, 198, 19, 175, 254, 202, 208, 0, 2, 0, 0,
    ];

    #[test]
    fn _count_kmer() {
        init();

        let file = generate_seq_file();

        let path = file.path().to_str().unwrap().to_string();

        let count = count_kmer(&[path], 5, 26).unwrap();

        let mut output = std::io::Cursor::new(Vec::new());
        count.serialize(&mut output).unwrap();

        assert_eq!(COUNT, output.into_inner());
    }

    #[test]
    fn _compute_abundance_threshold() {
        init();

        let output = vec![];
        let count = pcon::counter::Counter::deserialize(COUNT).unwrap();

        assert_eq!(2, compute_abundance_threshold(&count, output).unwrap());
    }

    static SOLID: &[u8] = &[
        31, 139, 8, 0, 0, 0, 0, 0, 4, 255, 77, 136, 193, 9, 0, 64, 8, 195, 94, 183, 255, 136, 55,
        134, 15, 145, 216, 10, 130, 208, 180, 52, 15, 40, 113, 147, 62, 223, 181, 248, 140, 93,
        161, 13, 1, 52, 40, 137, 69, 44, 65, 0, 0, 0,
    ];

    #[test]
    fn _read_solidity() {
        init();

        let count = pcon::counter::Counter::deserialize(COUNT).unwrap();
        let compute = pcon::solid::Solid::from_counter(&count, 2);

        let mut compute_out = vec![];
        compute.serialize(&mut compute_out).unwrap();

        assert_eq!(SOLID, compute_out.as_slice());

        let mut tmpfile = tempfile::NamedTempFile::new().unwrap();

        tmpfile.write_all(SOLID).unwrap();

        let path = tmpfile.path().to_str().unwrap().to_string();

        let read = read_solidity(path).unwrap();
        for kmer in 0..cocktail::kmer::get_kmer_space_size(5) {
            assert_eq!(compute.get(kmer), read.get(kmer));
        }
    }
}
