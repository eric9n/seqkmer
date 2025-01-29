use crate::fasta::{BufferFastaReader, FastaReader};
use crate::fastq::FastqReader;
use crate::reader::{detect_file_format, Reader};
use crate::seq::{Base, SeqFormat};
use crate::utils::OptionPair;
use std::io::Result;
use std::path::Path;

/// A reader for both FASTA and FASTQ files.
///
/// # Examples
///
/// ```
/// use seqkmer::{FastxReader, Reader, OptionPair};
/// use std::path::Path;
///
/// # fn main() -> std::io::Result<()> {
/// let path = Path::new("tests/data/test.fasta");
/// let mut reader = FastxReader::from_paths(OptionPair::Single(path), 0, 0)?;
///
/// while let Some(sequences) = reader.next()? {
///     for sequence in sequences {
///         println!("Sequence ID: {}", sequence.header.id);
///         println!("Sequence length: {}", sequence.body.single().unwrap().len());
///     }
/// }
/// # Ok(())
/// # }
/// ```
pub struct FastxReader<R: Reader> {
    inner: R,
}

impl<R: Reader> FastxReader<R> {
    /// Creates a new `FastxReader` with the given inner reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqkmer::{FastxReader, FastaReader};
    /// use std::path::Path;
    ///
    /// # fn main() -> std::io::Result<()> {
    /// let path = Path::new("tests/data/test.fasta");
    /// let fasta_reader = FastaReader::from_path(path, 0)?;
    /// let fastx_reader = FastxReader::new(fasta_reader);
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }
}

impl<R: Reader> Reader for FastxReader<R> {
    fn next(&mut self) -> Result<Option<Vec<Base<Vec<u8>>>>> {
        self.inner.next()
    }
}

impl FastxReader<Box<dyn Reader + Send>> {
    /// Creates a new `FastxReader` from file paths.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqkmer::{FastxReader, Reader, OptionPair};
    /// use std::path::Path;
    ///
    /// # fn main() -> std::io::Result<()> {
    /// let path = Path::new("tests/data/test.fasta");
    /// let reader = FastxReader::from_paths(OptionPair::Single(path), 0, 0)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_paths<P: AsRef<Path>>(
        paths: OptionPair<P>,
        file_index: usize,
        quality_score: i32,
    ) -> Result<Self> {
        let file_format = paths.map(|path: &P| detect_file_format(path));

        match file_format? {
            OptionPair::Single(SeqFormat::Fasta) => {
                let reader = FastaReader::from_path(paths.single().unwrap().as_ref(), file_index)?;
                Ok(Self::new(Box::new(reader) as Box<dyn Reader + Send>))
            }
            OptionPair::Single(SeqFormat::Fastq)
            | OptionPair::Pair(SeqFormat::Fastq, SeqFormat::Fastq) => {
                let reader = FastqReader::from_path(paths, file_index, quality_score)?;
                Ok(Self::new(Box::new(reader) as Box<dyn Reader + Send>))
            }
            _ => panic!("Unsupported file format combination"),
        }
    }

    /// Creates a new `FastxReader` using a buffered reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqkmer::{FastxReader, Reader, OptionPair};
    /// use std::path::Path;
    ///
    /// # fn main() -> std::io::Result<()> {
    /// let path = Path::new("tests/data/test.fasta");
    /// let reader = FastxReader::from_buffer_reader(OptionPair::Single(path), 0, 0)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_buffer_reader<P: AsRef<Path>>(
        paths: OptionPair<P>,
        file_index: usize,
        quality_score: i32,
    ) -> Result<Self> {
        let file_format = paths.map(|path: &P| detect_file_format(path));

        match file_format? {
            OptionPair::Single(SeqFormat::Fasta) => {
                let reader =
                    BufferFastaReader::from_path(paths.single().unwrap().as_ref(), file_index)?;
                Ok(Self::new(Box::new(reader) as Box<dyn Reader + Send>))
            }
            OptionPair::Single(SeqFormat::Fastq)
            | OptionPair::Pair(SeqFormat::Fastq, SeqFormat::Fastq) => {
                let reader = FastqReader::from_path(paths, file_index, quality_score)?;
                Ok(Self::new(Box::new(reader) as Box<dyn Reader + Send>))
            }
            _ => panic!("Unsupported file format combination"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn test_read_ont_fastq() -> std::io::Result<()> {
        let path = Path::new("tests/data/example_ont_reads.fastq");
        let mut reader = FastxReader::from_paths(OptionPair::Single(path), 0, 0)?;

        let mut read_count = 0;
        let mut total_bases = 0;

        while let Some(sequences) = reader.next()? {
            for sequence in sequences {
                read_count += 1;
                if let OptionPair::Single(seq) = sequence.body {
                    total_bases += seq.len();
                    
                    // Verify the first read's content
                    if read_count == 1 {
                        assert_eq!(sequence.header.id, "89a96608-1899-49e1-b077-767a40d5ae27");
                        
                        // Check sequence length
                        assert_eq!(seq.len(), 3928, "First sequence should be 3928 bases long");
                        
                        // Check sequence start and end
                        let start = std::str::from_utf8(&seq[..10]).unwrap();
                        let end = std::str::from_utf8(&seq[seq.len()-10..]).unwrap();
                        
                        assert_eq!(start, "ATGTTTTGTA");
                        assert_eq!(end, "GTGGTGCCAT");
                        
                        // Verify sequence only contains valid DNA characters
                        for &base in seq.iter() {
                            assert!(matches!(base, b'A' | b'T' | b'C' | b'G' | b'N'));
                        }
                    }
                }
            }
        }

        // Verify we read all 5 sequences from the file
        assert_eq!(read_count, 5, "Should have read 5 sequences");
        assert!(total_bases > 0, "Should have read some bases");
        
        // Verify total bases is reasonable
        assert!(total_bases > 10000, "Total bases should be substantial");

        Ok(())
    }


}
