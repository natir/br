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

use bv;
use std::io::Read;

fn get_k(path: &str) -> u8 {
    return (((std::fs::metadata(path).unwrap().len() as f64).log2() as u8 ) / 2) + 2;
}

pub fn read_existance(path: &str) -> (u8, bv::BitVec<u8>) {
    let mut reader = std::io::BufReader::new(std::fs::File::open(path).unwrap());

    let k = get_k(path);

    let mut data: Vec<u8> = vec![0u8; (pcon::io::read::get_kmer_space_size(k) / 8) as usize];
    reader.read_exact(&mut data).expect("Error durring reading of data");
    
    return (k, bv::BitVec::from_bits(&data));
}
