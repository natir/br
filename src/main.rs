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
use anyhow::Result;
use clap::Clap;

use br::*;

fn main() -> Result<()> {
    let mut params = cli::Command::parse();

    params = cli::check_params(params)?;

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

    let confirm = if let Some(val) = params.confirm {
        val
    } else {
        2
    };

    let max_search = if let Some(val) = params.max_search {
        val
    } else {
        7
    };

    if let Some(threads) = params.threads {
        log::info!("Set number of threads to {}", threads);

        set_nb_threads(threads);
    }

    let record_buffer = if let Some(len) = params.record_buffer {
        len
    } else {
        8192
    };

    let solid = cli::read_or_compute_solidity(
        params.solidity,
        params.kmer_solid,
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
