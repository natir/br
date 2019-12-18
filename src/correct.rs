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

#[derive(Debug, PartialEq)]
pub enum ErrorType {
    None = 0,
    Sub = 1,
    Del = 2,
    Ins = 4,
}


pub fn correct_read(seq: &[u8], k: u8, s: u8, data: &bv::BitVec<u8>) -> Vec<u8> {
    let mut corrected_seq: Vec<u8> = Vec::with_capacity(seq.len());
    
    if seq.len() < k as usize { return corrected_seq; } // if seq is lower than k we can't do anything

    for i in 0..(k-1) {
        corrected_seq.push(seq[i as usize]);
    }

    if seq.len() == k as usize { return corrected_seq; }

    let mut previous = false;
    for (i, mut kmer) in cocktail::tokenizer::Tokenizer::new(seq, k).map(|x| cocktail::kmer::cannonical(x, k)).enumerate() {
	debug!("position {}", i + k as usize);
	let nuc = seq[i + k as usize];
	
	let kmer_solidity = this_kmer_is_true(data, kmer, k);
        if !kmer_solidity && previous {
            debug!("kmer {} is not present", cocktail::kmer::kmer2seq(kmer, k));
	    if let Some((corrected_kmer, correct_seq)) = correct_kmer(data, kmer, k, s, seq, i + k as usize) {
                debug!("base write in correct read {:?}", correct_seq);
                kmer = corrected_kmer;
                for c_nuc in &correct_seq {
                    corrected_seq.push(*c_nuc);
                }
		previous = true;
            } else {
		previous = false;
                corrected_seq.push(nuc);
            }
        } else {
            // If kmer is present or if previous kmer isn't solid we didn't do anything
            debug!("kmer {} is present {} ", cocktail::kmer::kmer2seq(kmer, k), kmer_solidity);
	    previous = kmer_solidity;
            corrected_seq.push(nuc);
        }
    }

    return corrected_seq;
}

fn correct_kmer(data: &bv::BitVec<u8>, kmer: u64, k: u8, s: u8, seq: &[u8], i: usize) -> Option<(u64, Vec<u8>)> {
    let nuc = seq[i];
    let kmer_mask = (1 << k * 2) - 1;

    let mut scenario: Vec<(u8, u8, ErrorType, u64)> = Vec::new();

    let alt_nucs = get_alt_nucs(data, kmer, k);
    if alt_nucs.len() != 1 {
        debug!("multiple alternative kmer or no alternative {}", alt_nucs.len());
        return None;
    }
    
    for alt_nuc in alt_nucs {
        let alt_kmer = ((kmer >> 2) << 2) | (alt_nuc as u64);
	debug!("alt_kmer {}", cocktail::kmer::kmer2seq(alt_kmer, k));

        if let Some(limit) = get_end_of_subseq(i+1, s as usize, seq.len()) {
	    debug!("Test substitution");
            scenario.push((get_kmer_score(data, alt_kmer, k, &seq[i+1..limit], kmer_mask), alt_nuc, ErrorType::Sub, alt_kmer));
        }

        if let Some(limit) = get_end_of_subseq(i+2, s as usize, seq.len()) {
	    debug!("Test insertion");
            scenario.push((get_kmer_score(data, alt_kmer, k, &seq[i+2..limit], kmer_mask), alt_nuc, ErrorType::Ins, alt_kmer));
        }
        
        if let Some(limit) = get_end_of_subseq(i, s as usize, seq.len()) {
	    debug!("Test deletion");
            scenario.push((get_kmer_score(data, alt_kmer, k, &seq[i..limit], kmer_mask), alt_nuc, ErrorType::Del, alt_kmer));
        }
    }
    
    scenario.drain_filter(|x| x.0 != s);
    scenario.sort_by(|a, b| a.0.cmp(&b.0));
    
    if scenario.len() != 1 {
        debug!("No valid scenario or to many {:?}", scenario);
	return None;
    }
    
    let (score, alt_nuc, error_type, alt_kmer) = &scenario[0];
    
    return match error_type {
	ErrorType::Sub => Some((*alt_kmer, vec![cocktail::kmer::bit2nuc(*alt_nuc as u64)])),
        ErrorType::Ins => Some((alt_kmer >> 2, vec![])),
        ErrorType::Del => Some(((alt_kmer << 2) | cocktail::kmer::nuc2bit(nuc) as u64, vec![cocktail::kmer::bit2nuc(*alt_nuc as u64), nuc])),
        _ => None,
    };
}

fn get_end_of_subseq(begin: usize, offset: usize, max_length: usize) -> Option<usize> {
    return match begin + offset > max_length {
        false => Some(begin + offset),
        true  => None // we are at end of sequence
    };
}

fn get_kmer_score(data: &bv::BitVec<u8>, mut alt_kmer: u64, k: u8, nucs: &[u8], kmer_mask: u64) -> u8 {
    let mut score = 0;

    debug!("\t nucs to test {:?}", nucs);
    for nuc in nucs {
        alt_kmer = ((alt_kmer << 2) & kmer_mask) | cocktail::kmer::nuc2bit(*nuc);

	debug!("\t alt_kmer {} is present {}", cocktail::kmer::kmer2seq(alt_kmer, k), this_kmer_is_true(data, alt_kmer, k));
        if this_kmer_is_true(data, alt_kmer, k) {
            score += 1
        }
    }

    return score;
}


fn this_kmer_is_true(data: &bv::BitVec<u8>, kmer: u64, k: u8) -> bool {
    return data.get(cocktail::kmer::remove_first_bit(cocktail::kmer::cannonical(kmer, k)));
}

fn get_alt_nucs(data: &bv::BitVec<u8>, kmer: u64, k: u8) -> Vec<u8> {
    let mut correct_nuc: Vec<u8> = Vec::with_capacity(3);
    let kme = (kmer >> 2) << 2;

    for alt_nuc in (0..4).filter(|x| *x != (kmer & 0b11) as u8) {
        if this_kmer_is_true(data, kme ^ alt_nuc as u64, k) {
            correct_nuc.push(alt_nuc);
        }
    }

    return correct_nuc;
}

#[cfg(test)]
mod tests {
    
    use super::*;

    #[test]
    fn substitution() {
        let correct = b"ACTGACGA";
        let noisy   = b"ACTGATGA";

        let mut data: bv::BitVec<u8> = bv::BitVec::new_fill(false, 1 << 9);

        data.set(pcon::convert::hash(b"ACTGA", 5), true);
        data.set(pcon::convert::hash(b"CTGAC", 5), true);
        data.set(pcon::convert::hash(b"TGACG", 5), true);
        data.set(pcon::convert::hash(b"GACGA", 5), true);
        
        assert_eq!(correct, correct_read(noisy, 5, 2, &data).as_slice());
    }

    #[test]
    fn substitution_multiple_alternative_k() {
        let correct = b"ACTGACGA";
        let noisy   = b"ACTGATGA";

        let mut data: bv::BitVec<u8> = bv::BitVec::new_fill(false, 1 << 9);

        data.set(pcon::convert::hash(b"ACTGA", 5), true);
        data.set(pcon::convert::hash(b"CTGAC", 5), true);
        data.set(pcon::convert::hash(b"CTGAG", 5), true);
        data.set(pcon::convert::hash(b"TGACG", 5), true);
        data.set(pcon::convert::hash(b"GACGA", 5), true);
        
        assert_eq!(noisy, correct_read(noisy, 5, 2, &data).as_slice());
    }
    
    #[test]
    fn insertion() {
        let correct = b"ACTGACGA";
        let noisy   = b"ACTGATCGA";

        let mut data: bv::BitVec<u8> = bv::BitVec::new_fill(false, 1 << 9);

        data.set(pcon::convert::hash(b"ACTGA", 5), true);
        data.set(pcon::convert::hash(b"CTGAC", 5), true);
        data.set(pcon::convert::hash(b"TGACG", 5), true);
        data.set(pcon::convert::hash(b"GACGA", 5), true);
        
        assert_eq!(correct, correct_read(noisy, 5, 2, &data).as_slice());
    }

    #[test]
    fn insertion_multiple_alternative_k() {
        let correct = b"ACTGACGA";
        let noisy   = b"ACTGATCGA";

        let mut data: bv::BitVec<u8> = bv::BitVec::new_fill(false, 1 << 9);

        data.set(pcon::convert::hash(b"ACTGA", 5), true);
        data.set(pcon::convert::hash(b"CTGAC", 5), true);
        data.set(pcon::convert::hash(b"CTGAG", 5), true);
        data.set(pcon::convert::hash(b"TGACG", 5), true);
        data.set(pcon::convert::hash(b"GACGA", 5), true);
        
        assert_eq!(noisy, correct_read(noisy, 5, 2, &data).as_slice());
    }
    
    #[test]
    fn deletion() {
        let correct = b"ACTGACGA";
        let noisy   = b"ACTGAGA";

        let mut data: bv::BitVec<u8> = bv::BitVec::new_fill(false, 1 << 9);

        data.set(pcon::convert::hash(b"ACTGA", 5), true);
        data.set(pcon::convert::hash(b"CTGAC", 5), true);
        data.set(pcon::convert::hash(b"TGACG", 5), true);
        data.set(pcon::convert::hash(b"GACGA", 5), true);
        
        assert_eq!(correct, correct_read(noisy, 5, 2, &data).as_slice());
    }

    #[test]
    fn deletion_multiple_alternative_k() {
        let correct = b"ACTGACGA";
        let noisy   = b"ACTGAGA";

        let mut data: bv::BitVec<u8> = bv::BitVec::new_fill(false, 1 << 9);
        
        assert_eq!(noisy, correct_read(noisy, 5, 2, &data).as_slice());
    }
}
