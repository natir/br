//! Trait for hashset

/* std use */

/* crate use */

/* project use */

/* mod declaration */
pub mod hash;
pub mod pcon;

/* reexport */
pub use self::hash::Hash;
pub use self::pcon::Pcon;

pub trait KmerSet: Sync {
    fn get(&self, kmer: u64) -> bool;

    fn k(&self) -> u8;
}

pub type BoxKmerSet<'a> = Box<dyn KmerSet + 'a>;
