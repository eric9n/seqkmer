use crate::reader::{dyn_reader, Reader, BUFSIZE};
use crate::seq::{Base, SeqFormat, SeqHeader};
use crate::utils::OptionPair;
use std::io::{BufRead, BufReader, Read, Result};
use std::path::Path;

const SEQ_LIMIT: u64 = u64::pow(2, 32);

fn trim_line_end(buffer: &mut Vec<u8>) {
    while matches!(buffer.last(), Some(b'\n' | b'\r')) {
        buffer.pop();
    }
}

fn first_non_whitespace(bytes: &[u8]) -> Option<usize> {
    bytes.iter().position(|b| !b.is_ascii_whitespace())
}

fn is_comment(bytes: &[u8], start: usize) -> bool {
    matches!(bytes.get(start), Some(b';'))
}

fn append_seq_line(seq: &mut Vec<u8>, line: &[u8], start: usize) {
    for &base in &line[start..] {
        if !base.is_ascii_whitespace() {
            seq.push(base);
        }
    }
}

fn read_next_record<R: Read>(
    reader: &mut BufReader<R>,
    pending_header: &mut Option<Vec<u8>>,
    header: &mut Vec<u8>,
    seq: &mut Vec<u8>,
    line_buf: &mut Vec<u8>,
) -> Result<Option<usize>> {
    header.clear();
    seq.clear();

    if let Some(next_header) = pending_header.take() {
        header.extend_from_slice(&next_header);
    } else {
        loop {
            line_buf.clear();
            let bytes_read = reader.read_until(b'\n', line_buf)?;
            if bytes_read == 0 {
                return Ok(None);
            }
            trim_line_end(line_buf);
            let Some(start) = first_non_whitespace(line_buf) else {
                continue;
            };
            if is_comment(line_buf, start) {
                continue;
            }
            if line_buf[start] == b'>' {
                header.extend_from_slice(&line_buf[start..]);
                break;
            }
        }
    }

    let mut seq_lines = 0usize;

    loop {
        line_buf.clear();
        let bytes_read = reader.read_until(b'\n', line_buf)?;
        if bytes_read == 0 {
            break;
        }
        trim_line_end(line_buf);
        let Some(start) = first_non_whitespace(line_buf) else {
            continue;
        };
        if is_comment(line_buf, start) {
            continue;
        }
        if line_buf[start] == b'>' {
            pending_header.replace(line_buf[start..].to_vec());
            break;
        }
        append_seq_line(seq, line_buf, start);
        seq_lines += 1;
    }

    Ok(Some(seq_lines))
}

/// FastaReader for reading FASTA format files.
///
/// # Examples
///
/// ```
/// use seqkmer::{FastaReader, Reader};
/// use std::path::Path;
///
/// # fn main() -> std::io::Result<()> {
/// let path = Path::new("tests/data/test.fasta");
/// let mut reader = FastaReader::from_path(path, 0)?;
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
pub struct FastaReader<R>
where
    R: Read + Send,
{
    reader: BufReader<R>,
    file_index: usize,
    reads_index: usize,
    header: Vec<u8>,
    seq: Vec<u8>,
    line_buf: Vec<u8>,
    pending_header: Option<Vec<u8>>,

    // 批量读取
    batch_size: usize,
}

impl<R> FastaReader<R>
where
    R: Read + Send,
{
    /// Creates a new FastaReader with default capacity and batch size.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqkmer::FastaReader;
    /// use std::fs::File;
    ///
    /// # fn main() -> std::io::Result<()> {
    /// let file = File::open("tests/data/test.fasta")?;
    /// let reader = FastaReader::new(file, 0);
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(reader: R, file_index: usize) -> Self {
        Self::with_capacity(reader, file_index, BUFSIZE, 30)
    }

    /// Creates a new FastaReader with specified capacity and batch size.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqkmer::FastaReader;
    /// use std::fs::File;
    ///
    /// # fn main() -> std::io::Result<()> {
    /// let file = File::open("tests/data/test.fasta")?;
    /// let reader = FastaReader::with_capacity(file, 0, 4096, 50);
    /// # Ok(())
    /// # }
    /// ```
    pub fn with_capacity(reader: R, file_index: usize, capacity: usize, batch_size: usize) -> Self {
        assert!(capacity >= 3);
        Self {
            reader: BufReader::with_capacity(capacity, reader),
            file_index,
            reads_index: 0,
            header: Vec::new(),
            seq: Vec::new(),
            line_buf: Vec::new(),
            pending_header: None,
            batch_size,
        }
    }

    pub fn read_next(&mut self) -> Result<Option<()>> {
        read_next_record(
            &mut self.reader,
            &mut self.pending_header,
            &mut self.header,
            &mut self.seq,
            &mut self.line_buf,
        )
        .map(|opt| opt.map(|_| ()))
    }

    pub fn _next(&mut self) -> Result<Option<(usize, Base<Vec<u8>>)>> {
        if self.read_next()?.is_none() {
            return Ok(None);
        }

        let seq_len = self.seq.len();
        // 检查seq的长度是否大于2的32次方
        if seq_len as u64 > SEQ_LIMIT {
            eprintln!("Sequence length exceeds 2^32, which is not handled.");
            return Ok(None);
        }

        let seq_id = unsafe {
            let slice = if self.header.starts_with(b">") {
                &self.header[1..]
            } else {
                &self.header[..]
            };

            let s = std::str::from_utf8_unchecked(slice);
            let first_space_index = s
                .find(|c: char| c.is_whitespace() || c == '\u{1}')
                .unwrap_or(s.len());

            // 直接从原始切片创建第一个单词的切片
            &s[..first_space_index]
        };
        self.reads_index += 1;

        let seq_header = SeqHeader {
            file_index: self.file_index,
            reads_index: self.reads_index,
            format: SeqFormat::Fasta,
            id: seq_id.to_owned(),
        };
        Ok(Some((
            seq_len,
            Base::new(seq_header, OptionPair::Single(self.seq.to_owned())),
        )))
    }
}

impl FastaReader<Box<dyn Read + Send>> {
    /// Creates a new FastaReader from a file path.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqkmer::FastaReader;
    /// use std::path::Path;
    ///
    /// # fn main() -> std::io::Result<()> {
    /// let path = Path::new("tests/data/test.fasta");
    /// let reader = FastaReader::from_path(path, 0)?;
    /// # Ok(())
    /// # }
    /// ```
    #[inline]
    pub fn from_path<P: AsRef<Path>>(path: P, file_index: usize) -> Result<Self> {
        let reader = dyn_reader(path)?;
        Ok(Self::new(reader, file_index))
    }
}

impl<R: Read + Send> Reader for FastaReader<R> {
    fn next(&mut self) -> Result<Option<Vec<Base<Vec<u8>>>>> {
        let mut seqs = Vec::new();
        let mut total_bytes = 0;
        let max_bytes = 10 * 1024 * 1024;

        for _ in 0..self.batch_size {
            if let Some((seq_len, seq)) = self._next()? {
                seqs.push(seq);
                total_bytes += seq_len;
                if total_bytes > max_bytes {
                    break;
                }
            } else {
                break;
            }
        }

        Ok(if seqs.is_empty() { None } else { Some(seqs) })
    }
}

/// BufferFastaReader for reading FASTA format files with buffering.
///
/// # Examples
///
/// ```
/// use seqkmer::{BufferFastaReader, Reader};
/// use std::path::Path;
///
/// # fn main() -> std::io::Result<()> {
/// let path = Path::new("tests/data/test.fasta");
/// let mut reader = BufferFastaReader::from_path(path, 0)?;
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
pub struct BufferFastaReader<R>
where
    R: Read + Send,
{
    reader: BufReader<R>,
    file_index: usize,
    reads_index: usize,
    header: Vec<u8>,
    seq: Vec<u8>,
    line_buf: Vec<u8>,
    pending_header: Option<Vec<u8>>,
    line_num: usize,

    // 批量读取
    batch_size: usize,
}

impl<R> BufferFastaReader<R>
where
    R: Read + Send,
{
    /// Creates a new BufferFastaReader with default capacity and batch size.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqkmer::BufferFastaReader;
    /// use std::fs::File;
    ///
    /// # fn main() -> std::io::Result<()> {
    /// let file = File::open("tests/data/test.fasta")?;
    /// let reader = BufferFastaReader::new(file, 0);
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(reader: R, file_index: usize) -> Self {
        Self::with_capacity(reader, file_index, BUFSIZE, 60)
    }

    /// Creates a new BufferFastaReader with specified capacity and batch size.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqkmer::BufferFastaReader;
    /// use std::fs::File;
    ///
    /// # fn main() -> std::io::Result<()> {
    /// let file = File::open("tests/data/test.fasta")?;
    /// let reader = BufferFastaReader::with_capacity(file, 0, 4096, 50);
    /// # Ok(())
    /// # }
    /// ```
    pub fn with_capacity(reader: R, file_index: usize, capacity: usize, batch_size: usize) -> Self {
        assert!(capacity >= 3);
        Self {
            reader: BufReader::with_capacity(capacity, reader),
            file_index,
            reads_index: 0,
            header: Vec::new(),
            seq: Vec::new(),
            line_buf: Vec::new(),
            pending_header: None,
            line_num: 0,
            batch_size,
        }
    }

    pub fn read_next(&mut self) -> Result<Option<()>> {
        self.line_num = 0;
        read_next_record(
            &mut self.reader,
            &mut self.pending_header,
            &mut self.header,
            &mut self.seq,
            &mut self.line_buf,
        )
        .map(|opt| opt.map(|count| {
            self.line_num = count;
        }))
    }

    pub fn _next(&mut self) -> Result<Option<Base<Vec<u8>>>> {
        if self.read_next()?.is_none() {
            return Ok(None);
        }

        let seq_len = self.seq.len();
        // 检查seq的长度是否大于2的32次方
        if seq_len as u64 > SEQ_LIMIT {
            eprintln!("Sequence length exceeds 2^32, which is not handled.");
            return Ok(None);
        }

        let seq_id = unsafe {
            let slice = if self.header.starts_with(b">") {
                &self.header[1..]
            } else {
                &self.header[..]
            };

            let s = std::str::from_utf8_unchecked(slice);
            let first_space_index = s
                .find(|c: char| c.is_whitespace() || c == '\u{1}')
                .unwrap_or(s.len());
            // let first_space_index = s
            //     .as_bytes()
            //     .iter()
            //     .position(|&c| c == b' ')
            //     .unwrap_or(s.len());

            // 直接从原始切片创建第一个单词的切片
            &s[..first_space_index]
        };
        self.reads_index += 1;

        let seq_header = SeqHeader {
            file_index: self.file_index,
            reads_index: self.reads_index,
            format: SeqFormat::Fasta,
            id: seq_id.to_owned(),
        };
        Ok(Some(Base::new(
            seq_header,
            OptionPair::Single(self.seq.to_owned()),
        )))
    }
}

impl BufferFastaReader<Box<dyn Read + Send>> {
    /// Creates a new BufferFastaReader from a file path.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqkmer::BufferFastaReader;
    /// use std::path::Path;
    ///
    /// # fn main() -> std::io::Result<()> {
    /// let path = Path::new("tests/data/test.fasta");
    /// let reader = BufferFastaReader::from_path(path, 0)?;
    /// # Ok(())
    /// # }
    /// ```
    #[inline]
    pub fn from_path<P: AsRef<Path>>(path: P, file_index: usize) -> Result<Self> {
        let reader = dyn_reader(path)?;
        Ok(Self::new(reader, file_index))
    }
}

impl<R: Read + Send> Reader for BufferFastaReader<R> {
    fn next(&mut self) -> Result<Option<Vec<Base<Vec<u8>>>>> {
        let mut seqs = Vec::new();
        for _ in 0..self.batch_size {
            match self._next()? {
                Some(seq) => seqs.push(seq),
                None => break,
            }
        }

        Ok(if seqs.is_empty() { None } else { Some(seqs) })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn collect_all<R: Read + Send, F: FnOnce(R) -> T, T: Reader>(
        reader_builder: F,
        source: R,
    ) -> Result<Vec<Base<Vec<u8>>>> {
        let mut reader = reader_builder(source);
        let mut collected = Vec::new();
        while let Some(batch) = reader.next()? {
            collected.extend(batch);
        }
        Ok(collected)
    }

    #[test]
    fn fasta_reader_handles_multiline_sequences() -> Result<()> {
        let data = b">seq1 desc\nACGT\nA T G C \n; comment line\n>seq2\nNNNN\r\n\r\nNN\tNN\n>seq3\n"
            .to_vec();
        let cursor = Cursor::new(data.clone());
        let mut reader = FastaReader::with_capacity(cursor, 0, 128, 4);
        let mut collected = Vec::new();
        while let Some(batch) = reader.next()? {
            collected.extend(batch);
        }

        assert_eq!(collected.len(), 3);

        let seq1 = collected[0].body.single().unwrap();
        assert_eq!(collected[0].header.id, "seq1");
        assert_eq!(seq1.as_slice(), b"ACGTATGC");

        let seq2 = collected[1].body.single().unwrap();
        assert_eq!(collected[1].header.id, "seq2");
        assert_eq!(seq2.as_slice(), b"NNNNNNNN");

        let seq3 = collected[2].body.single().unwrap();
        assert_eq!(collected[2].header.id, "seq3");
        assert!(seq3.is_empty());

        Ok(())
    }

    #[test]
    fn buffer_fasta_reader_matches_standard_reader() -> Result<()> {
        let data = b">seq1 desc\nACGT\nA T G C \n; comment line\n>seq2\nNNNN\r\n\r\nNN\tNN\n>seq3\n"
            .to_vec();

        let collected_fast = collect_all(
            |reader| FastaReader::with_capacity(reader, 0, 128, 4),
            Cursor::new(data.clone()),
        )?;
        let collected_buffered = collect_all(
            |reader| BufferFastaReader::with_capacity(reader, 0, 128, 4),
            Cursor::new(data),
        )?;

        assert_eq!(collected_fast.len(), collected_buffered.len());
        for (a, b) in collected_fast.iter().zip(collected_buffered.iter()) {
            assert_eq!(a.header.id, b.header.id);
            assert_eq!(a.body.single().unwrap(), b.body.single().unwrap());
        }

        Ok(())
    }
}
