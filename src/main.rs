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

#![feature(drain_filter)]

/* crate declaration */
extern crate bio;
extern crate clap;
extern crate bv;
extern crate cocktail;
#[macro_use] extern crate log;

/* crate use */
use clap::{App, Arg};

/* local mod */
mod correct;

fn main() {
    env_logger::builder()
	.format_timestamp(None)
	.init();
    
    let matches = App::new("br")
        .version("0.1")
        .author("Pierre Marijon <pmarijon@mmci.uni-saarland.de>")
        .about("A simple long read correcteur based on kmer spectrum method")
        .arg(Arg::with_name("solidity")
             .short("s")
             .long("solidity")
             .required(true)
             .takes_value(true)
             .help("solidity bitfield produce by Pcon")
        )
        .arg(Arg::with_name("input")
             .short("i")
             .long("input")
             .required(true)
             .takes_value(true)
             .multiple(true)
             .help("fasta file to be correct")
        )
        .arg(Arg::with_name("output")
             .short("o")
             .long("output")
             .required(true)
             .takes_value(true)
             .multiple(true)
             .help("path where corrected read was write")
        )
        .arg(Arg::with_name("confirm")
             .short("c")
             .long("confirm")
             .default_value("2")
             .takes_value(true)
             .help("number of kmer need to be solid after correction to validate it")
        )
        .get_matches();


    let inputs_path: Vec<&str> = matches.values_of("input").unwrap().collect();
    let outputs_path: Vec<&str> = matches.values_of("output").unwrap().collect();
    let solidity = matches.value_of("solidity").unwrap().parse::<u8>().expect("We can't parse the parameter solidity");

    let solidity_path = matches.value_of("solidity").unwrap();    
    let solidity_file = std::fs::File::open(solidity_path).expect("Error when we try to open solidity file");
    let file_size = solidity_file.metadata().expect("Error we can't access to metadat").len();
    let (k, data) = cocktail::io::read_solidity_bitfield(std::io::BufReader::new(solidity_file), file_size);

    if inputs_path.len() != outputs_path.len() {
        eprintln!("number of inputs isn't equal to number of outputs.");
        return ();
    }

    for (input_path, output_path) in inputs_path.iter().zip(outputs_path) {
        let reader = bio::io::fasta::Reader::new(std::io::BufReader::new(std::fs::File::open(input_path).unwrap()));
        let mut write = bio::io::fasta::Writer::new(std::io::BufWriter::new(std::fs::File::create(output_path).unwrap()));

        for record in reader.records() {
            let result = record.unwrap();
            let seq = result.seq();
            
            let mut right_left = correct::correct_read(seq, k, solidity, &data);

            //right_left.reverse();
            //let mut left_right = correct::correct_read(&right_left, k, solidity, &data);
            //left_right.reverse();
            
            let left_right = right_left;
            
            write.write_record(&bio::io::fasta::Record::with_attrs(result.id(), result.desc(), &left_right)).expect("Error when we try to write corrected sequence");
        }
    }
}
