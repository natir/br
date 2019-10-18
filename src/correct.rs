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

enum_from_primitive!{
    #[derive(Debug, PartialEq)]
    pub enum ErrorType {
        None = 0,
        Subs = 1,
        Dele = 2,
        SubsDele = 3,
        Inse = 4,
        SubsInse = 5,
        DeleInse = 6,
        SubsDeleInse = 7,
    }
}

pub fn correct_read(seq: &[u8], k: u8, s: u8, data: &bv::BitVec<u8>) -> Vec<u8> {
    let mut corrected_seq: Vec<u8> = Vec::with_capacity(seq.len());
    
    if seq.len() < k as usize { return corrected_seq; } // if seq is lower than k we can't do anything

    let kmer_mask = (1 << k * 2) - 1;

    for i in 0..(k-1) {
        corrected_seq.push(seq[i as usize]);
    }

    let mut prev_is_solid = true;
    
    let mut kmer = pcon::convert::seq2bit(&seq[0..(k-1) as usize]);

    println!("kmer {:?} bit {:010b}", &seq[0..(k-1) as usize], kmer);
    
    if seq.len() == k as usize { return corrected_seq; }
    
    for i in ((k - 1) as usize)..(seq.len()) {

        let nuc = seq[i];
        kmer = ((kmer << 2) & kmer_mask) | pcon::convert::nuc2bit(nuc);
        print!(" kmer {:010b}", kmer);
        if !this_kmer_is_true(data, kmer, k) && prev_is_solid {
            if let Some(kmer) = correct_kmer(data, kmer, k, s, seq, &mut corrected_seq, i) {
                prev_is_solid = true;
            } else {
                prev_is_solid = false;
            }
        } else {
            // If kmer is present we didn't do anything
            print!(" not error");
            corrected_seq.push(nuc);
            prev_is_solid = true;
        }

        println!("");
    }

    return corrected_seq;
}

fn correct_kmer(data: &bv::BitVec<u8>, mut kmer: u64, k: u8, s: u8, seq: &[u8], corrected_seq: &mut Vec<u8>, i: usize) -> Option<u64> {
    let mut error_type = 1;
    let nuc = seq[i];
    let kmer_mask = (1 << k * 2) - 1;

    print!(" {}", nuc as char);

    // This kmer didn't exist in solid kmer
    print!(" error");
    if let Some(alt_nuc) = get_good_kmer(data, kmer, k) {
        // If we replace last nuc by another one we found a solid kmer
        print!(" alt_nuc {:02b}", alt_nuc);
        let alt_kmer = ((kmer >> 2) << 2) | (alt_nuc as u64);
        print!(" alt_kmer {:010b}", alt_kmer);                

        // We try to see if by replace the nuc by the new one in next s kmer, there kmer was solid two 
        if let Some(limit) = get_end_of_subseq(i+1, s as usize, seq.len()) {
            // Substitution
            if next_kmer_was_true(data, alt_kmer, k, &seq[(i+1)..limit], kmer_mask) {
                error_type += 2;
            }
        }

        // If the next nuc was equal to the alt_nuc this is a deletion we try to validate s next kmer with this correction
        // TODO check
        if let Some(limit) = get_end_of_subseq(i+2, s as usize, seq.len()) {
            // Insertition
            if alt_nuc as u64 == pcon::convert::nuc2bit(seq[i+1]) {
                if next_kmer_was_true(data, alt_kmer, k, &seq[(i+2)..limit], kmer_mask) {
                    error_type += 4;
                }
            }
        }

        // We try to see if by consider this error was an insertion the next s kmer, there kmer was solid two 
        if let Some(limit) = get_end_of_subseq(i, s as usize, seq.len()) {
            // Deletion
            if next_kmer_was_true(data, alt_kmer, k, &seq[(i)..limit], kmer_mask) {
                error_type += 8;
            }
        }

        // If it's only a Substitution a Deletion or a Insertion we perform correction
        print!(" error type {}", error_type);
        print!(" kmer {:010b}", kmer);
        match error_type {
            3 => {
                corrected_seq.push(pcon::convert::bit2nuc(alt_nuc as u64));
                kmer = alt_kmer;
            }
            5 => {
                kmer = alt_kmer >> 2;
            }
            9 => {
                corrected_seq.push(pcon::convert::bit2nuc(alt_nuc as u64));
                corrected_seq.push(pcon::convert::bit2nuc(nuc as u64));
                kmer = (alt_kmer << 2) | nuc as u64;
            },
            _ => (),
        }
        print!(" new kmer {:010b}", kmer);

        return Some(kmer);
    }
    else {
        // If we have some possiblity to correcte kmer we didn't do anything
        corrected_seq.push(nuc);

        return None;
    }
}



fn get_end_of_subseq(begin: usize, offset: usize, max_length: usize) -> Option<usize> {
    println!("\n\tbegin {} offset {} begin + offset {} max_length {}", begin, offset, begin + offset, max_length);
    return match begin + offset > max_length {
        false => Some(begin + offset),
        true  => None // we are at end of sequence
    };
}

fn next_kmer_was_true(data: &bv::BitVec<u8>, mut alt_kmer: u64, k: u8, nucs: &[u8], kmer_mask: u64) -> bool {
    print!(" next nucs {:?}", nucs);
    for nuc in nucs {
        alt_kmer = ((alt_kmer << 2) & kmer_mask) | pcon::convert::nuc2bit(*nuc);

        if !this_kmer_is_true(data, alt_kmer, k) {
            return false;
        }
    }

    return true;
}


fn this_kmer_is_true(data: &bv::BitVec<u8>, kmer: u64, k: u8) -> bool {
    return data.get(pcon::convert::remove_first_bit(pcon::convert::cannonical(kmer, k)));
}

fn get_good_kmer(data: &bv::BitVec<u8>, kmer: u64, k: u8) -> Option<u8> {
    let mut correct_kmer: Vec<u8> = Vec::with_capacity(4);
    let kme = (kmer >> 2) << 2;

    println!("");
    for alt_nuc in (0..4).filter(|x| *x != (kmer & 0b11) as u8) {
        println!("\tpossible alt_nuc {}", alt_nuc);
        println!("\tsubk {:010b} new kmer {:010b}", kme, kme ^ alt_nuc as u64);
        if this_kmer_is_true(data, kme ^ alt_nuc as u64, k) {
            correct_kmer.push(alt_nuc);
        }
    }
    println!("\t{:?}", correct_kmer);
    return match correct_kmer.len() {
        1 => Some(correct_kmer.pop().expect("Impossible error 1 contact author")),
        _ => None,
    };
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

        println!("ACTGA val {:010b}", pcon::convert::seq2bit(b"ACTGA"));
        println!("CTGAC val {:010b}", pcon::convert::seq2bit(b"CTGAC"));
        println!("CTGAT val {:010b}", pcon::convert::seq2bit(b"CTGAT"));
        println!("TGACG val {:010b}", pcon::convert::seq2bit(b"TGACG"));
        println!("GACGA val {:010b}", pcon::convert::seq2bit(b"GACGA"));


        println!("ACTGA exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"ACTGA"), 5));
        println!("CTGAC exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"CTGAC"), 5));
        println!("TGACG exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"TGACG"), 5));
        println!("GACGA exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"GACGA"), 5));
        println!("CTGAT exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"CTGAT"), 5));
        
        
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

        println!("ACTGA val {:010b}", pcon::convert::seq2bit(b"ACTGA"));
        println!("CTGAC val {:010b}", pcon::convert::seq2bit(b"CTGAC"));
        println!("CTGAT val {:010b}", pcon::convert::seq2bit(b"CTGAT"));
        println!("TGACG val {:010b}", pcon::convert::seq2bit(b"TGACG"));
        println!("GACGA val {:010b}", pcon::convert::seq2bit(b"GACGA"));


        println!("ACTGA exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"ACTGA"), 5));
        println!("CTGAC exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"CTGAC"), 5));
        println!("TGACG exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"TGACG"), 5));
        println!("GACGA exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"GACGA"), 5));
        println!("CTGAT exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"CTGAT"), 5));
        
        
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

        println!("ACTGA val {:010b}", pcon::convert::seq2bit(b"ACTGA"));
        println!("CTGAC val {:010b}", pcon::convert::seq2bit(b"CTGAC"));
        println!("CTGAT val {:010b}", pcon::convert::seq2bit(b"CTGAT"));
        println!("TGACG val {:010b}", pcon::convert::seq2bit(b"TGACG"));
        println!("GACGA val {:010b}", pcon::convert::seq2bit(b"GACGA"));


        println!("ACTGA exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"ACTGA"), 5));
        println!("CTGAC exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"CTGAC"), 5));
        println!("TGACG exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"TGACG"), 5));
        println!("GACGA exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"GACGA"), 5));
        println!("CTGAT exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"CTGAT"), 5));
        
        
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

        println!("ACTGA val {:010b}", pcon::convert::seq2bit(b"ACTGA"));
        println!("CTGAC val {:010b}", pcon::convert::seq2bit(b"CTGAC"));
        println!("CTGAT val {:010b}", pcon::convert::seq2bit(b"CTGAT"));
        println!("TGACG val {:010b}", pcon::convert::seq2bit(b"TGACG"));
        println!("GACGA val {:010b}", pcon::convert::seq2bit(b"GACGA"));


        println!("ACTGA exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"ACTGA"), 5));
        println!("CTGAC exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"CTGAC"), 5));
        println!("TGACG exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"TGACG"), 5));
        println!("GACGA exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"GACGA"), 5));
        println!("CTGAT exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"CTGAT"), 5));
        
        
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

        println!("ACTGA val {:010b}", pcon::convert::seq2bit(b"ACTGA"));
        println!("CTGAC val {:010b}", pcon::convert::seq2bit(b"CTGAC"));
        println!("CTGAG val {:010b}", pcon::convert::seq2bit(b"CTGAG"));
        println!("TGACG val {:010b}", pcon::convert::seq2bit(b"TGACG"));
        println!("GACGA val {:010b}", pcon::convert::seq2bit(b"GACGA"));


        println!("ACTGA exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"ACTGA"), 5));
        println!("CTGAC exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"CTGAC"), 5));
        println!("TGACG exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"TGACG"), 5));
        println!("GACGA exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"GACGA"), 5));
        println!("CTGAT exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"CTGAT"), 5));
        
        
        assert_eq!(correct, correct_read(noisy, 5, 2, &data).as_slice());
    }

    #[test]
    fn deletion_multiple_alternative_k() {
        let correct = b"ACTGACGA";
        let noisy   = b"ACTGAGA";

        let mut data: bv::BitVec<u8> = bv::BitVec::new_fill(false, 1 << 9);

        data.set(pcon::convert::hash(b"ACTGA", 5), true);
        data.set(pcon::convert::hash(b"CTGAC", 5), true);
        data.set(pcon::convert::hash(b"CTGAT", 5), true);
        data.set(pcon::convert::hash(b"TGACG", 5), true);
        data.set(pcon::convert::hash(b"GACGA", 5), true);

        println!("ACTGA val {:010b}", pcon::convert::seq2bit(b"ACTGA"));
        println!("CTGAC val {:010b}", pcon::convert::seq2bit(b"CTGAC"));
        println!("CTGAG val {:010b}", pcon::convert::seq2bit(b"CTGAG"));
        println!("TGACG val {:010b}", pcon::convert::seq2bit(b"TGACG"));
        println!("GACGA val {:010b}", pcon::convert::seq2bit(b"GACGA"));


        println!("ACTGA exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"ACTGA"), 5));
        println!("CTGAC exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"CTGAC"), 5));
        println!("TGACG exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"TGACG"), 5));
        println!("GACGA exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"GACGA"), 5));
        println!("CTGAT exist {}", this_kmer_is_true(&data, pcon::convert::seq2bit(b"CTGAT"), 5));
        
        
        assert_eq!(noisy, correct_read(noisy, 5, 2, &data).as_slice());
    }
}
