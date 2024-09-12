use crate::utils::OptionPair;

/// Represents the format of a sequence file.
///
/// # Examples
///
/// ```
/// use seqkmer::SeqFormat;
///
/// let format = SeqFormat::Fasta;
/// assert_eq!(format, SeqFormat::Fasta);
///
/// let format = SeqFormat::Fastq;
/// assert_eq!(format, SeqFormat::Fastq);
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub enum SeqFormat {
    Fasta,
    Fastq,
}

/// Represents the header information of a sequence.
///
/// # Examples
///
/// ```
/// use seqkmer::{SeqHeader, SeqFormat};
///
/// let header = SeqHeader {
///     id: "seq1".to_string(),
///     file_index: 0,
///     reads_index: 1,
///     format: SeqFormat::Fasta,
/// };
///
/// assert_eq!(header.id, "seq1");
/// assert_eq!(header.file_index, 0);
/// assert_eq!(header.reads_index, 1);
/// assert_eq!(header.format, SeqFormat::Fasta);
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeqHeader {
    pub id: String,
    pub file_index: usize,
    pub reads_index: usize,
    pub format: SeqFormat,
}

/// Represents a base structure containing a header and a body.
///
/// # Examples
///
/// ```
/// use seqkmer::{Base, SeqHeader, SeqFormat, OptionPair};
///
/// let header = SeqHeader {
///     id: "seq1".to_string(),
///     file_index: 0,
///     reads_index: 1,
///     format: SeqFormat::Fasta,
/// };
///
/// let body = OptionPair::Single(vec![65, 84, 67, 71]); // "ATCG"
///
/// let base = Base::new(header, body);
///
/// assert_eq!(base.header.id, "seq1");
/// assert_eq!(base.body.single().unwrap(), &vec![65, 84, 67, 71]);
/// ```
#[derive(Debug)]
pub struct Base<T> {
    pub header: SeqHeader,
    pub body: OptionPair<T>,
}

impl<T> Base<T> {
    /// Creates a new Base instance.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqkmer::{Base, SeqHeader, SeqFormat, OptionPair};
    ///
    /// let header = SeqHeader {
    ///     id: "seq1".to_string(),
    ///     file_index: 0,
    ///     reads_index: 1,
    ///     format: SeqFormat::Fasta,
    /// };
    ///
    /// let body = OptionPair::Single(vec![65, 84, 67, 71]); // "ATCG"
    ///
    /// let base = Base::new(header, body);
    ///
    /// assert_eq!(base.header.id, "seq1");
    /// assert_eq!(base.body.single().unwrap(), &vec![65, 84, 67, 71]);
    /// ```
    pub fn new(header: SeqHeader, body: OptionPair<T>) -> Self {
        Self { header, body }
    }

    /// Maps the body of the Base instance using a provided function.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqkmer::{Base, SeqHeader, SeqFormat, OptionPair};
    ///
    /// let header = SeqHeader {
    ///     id: "seq1".to_string(),
    ///     file_index: 0,
    ///     reads_index: 1,
    ///     format: SeqFormat::Fasta,
    /// };
    ///
    /// let body = OptionPair::Single(vec![65, 84, 67, 71]); // "ATCG"
    ///
    /// let base = Base::new(header, body);
    ///
    /// let mapped_base = base.map(|v| Ok::<_, ()>(v.len())).unwrap();
    ///
    /// assert_eq!(mapped_base.header.id, "seq1");
    /// assert_eq!(mapped_base.body.single().unwrap(), &4);
    /// ```
    pub fn map<U, E, F>(&self, mut f: F) -> Result<Base<U>, E>
    where
        F: FnMut(&T) -> Result<U, E>,
    {
        self.body.map(|t| f(&t)).map(|body| Base {
            header: self.header.clone(),
            body,
        })
    }
}
