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
