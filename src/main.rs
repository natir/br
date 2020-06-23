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
use structopt::StructOpt;
use anyhow::{anyhow, Result, Context};

use br::*;
use br::error::*;
use br::error::IO::*;

fn main() -> Result<()> {
    env_logger::builder()
	.format_timestamp(None)
	.init();
    
    let mut params = cli::Command::from_args();

    params = cli::check_params(params)?;

    let solidity_reader = std::io::BufReader::new(
	std::fs::File::open(&params.solidity)
	    .with_context(|| {
		Error::IO(CantOpenFile)
	    })
	    .with_context(|| {
		anyhow!("File {}", params.solidity.clone())
	    })?
    );
    let solid = pcon::solid::Solid::deserialize(solidity_reader)?;
    let k = solid.k;

    correct::init_masks(k);
    
    for (input, output) in params.inputs.iter().zip(params.outputs) {
        let reader = bio::io::fasta::Reader::new(
	    std::io::BufReader::new(
		std::fs::File::open(input)
		    .with_context(|| {
			Error::IO(CantOpenFile)
		    })
		    .with_context(|| {
			anyhow!("File {}", input.clone())
		    })?
	    )
	);
	
        let mut write = bio::io::fasta::Writer::new(
	    std::io::BufWriter::new(
		std::fs::File::create(&output)
		    .with_context(|| {
			Error::IO(CantCreateFile)
		    })
		    .with_context(|| {
			anyhow!("File {}", output.clone())
		    })?
	    )
	);

        for record in reader.records() {
            let result = record.unwrap();
            let seq = result.seq();
            
            let correct = correct::correct_seq(seq, params.confirm, &solid);

            write.write_record(&bio::io::fasta::Record::with_attrs(result.id(), result.desc(), &correct))
		.with_context(|| {
		    Error::IO(IO::ErrorDurringWrite)
		})
		.with_context(|| {
		    anyhow!("File {}", input.clone())
		})?;
        }
    }

    Ok(())
}
