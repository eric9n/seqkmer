// Modules and public exports
pub mod fasta;
pub mod fastq;
pub mod fastx;
pub mod feat;
pub mod mmscanner;
pub mod parallel;
pub mod reader;
pub mod seq;
pub mod utils;

pub use fasta::BufferFastaReader;
pub use fasta::FastaReader;
pub use fastq::FastqReader;
pub use fastx::FastxReader;
pub use feat::constants::*;
pub use feat::*;
pub use mmscanner::{scan_sequence, Cursor, MinimizerData, MinimizerIterator, MinimizerWindow};
pub use parallel::create_reader;
pub use parallel::{
    buffer_map_parallel, buffer_read_parallel, read_parallel, ParallelItem, ParallelResult,
};
pub use reader::*;
pub use seq::{Base, SeqFormat, SeqHeader};
pub use utils::OptionPair;
