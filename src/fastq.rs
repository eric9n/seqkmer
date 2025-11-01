use crate::reader::{dyn_reader, trim_end, trim_pair_info, Reader, BUFSIZE};
use crate::seq::{Base, SeqFormat, SeqHeader};
use crate::utils::OptionPair;
use std::io::{BufRead, BufReader, Read, Result};
use std::path::Path;


/// 读取模式，用于判断是单端读取还是双端读取
enum SingleReadMode {
    Unknown,
    Interleaved,
    SingleEnd,
}


struct QReader<R: Read + Send> {
    reader: BufReader<R>,
    quality_score: i32,

    header: Vec<u8>,
    seq: Vec<u8>,
    plus: Vec<u8>,
    quals: Vec<u8>,
}

impl<R> QReader<R>
where
    R: Read + Send,
{
    pub fn with_capacity(reader: R, capacity: usize, quality_score: i32) -> Self {
        assert!(capacity >= 3);
        Self {
            reader: BufReader::with_capacity(capacity, reader),
            header: Vec::new(),
            seq: Vec::new(),
            plus: Vec::new(),
            quals: Vec::new(),
            quality_score,
        }
    }

    pub fn read_next(&mut self) -> Result<Option<()>> {
        // 读取fastq文件header部分
        self.header.clear();
        if self.reader.read_until(b'\n', &mut self.header)? == 0 {
            return Ok(None);
        }
        // 读取fastq文件seq部分
        self.seq.clear();
        if self.reader.read_until(b'\n', &mut self.seq)? == 0 {
            return Ok(None);
        }
        trim_end(&mut self.seq);

        // 读取fastq文件+部分
        self.plus.clear();
        if self.reader.read_until(b'\n', &mut self.plus)? == 0 {
            return Ok(None);
        }

        // 读取fastq文件quals部分
        self.quals.clear();
        if self.reader.read_until(b'\n', &mut self.quals)? == 0 {
            return Ok(None);
        }
        trim_end(&mut self.quals);

        if self.quality_score > 0 {
            for (base, &qscore) in self.seq.iter_mut().zip(self.quals.iter()) {
                if (qscore as i32 - '!' as i32) < self.quality_score {
                    *base = b'x';
                }
            }
        }

        Ok(Some(()))
    }
}

/// FastqReader for reading FASTQ format files.
///
/// # Examples
///
/// ```
/// use seqkmer::{FastqReader, Reader, OptionPair};
/// use std::path::Path;
///
/// # fn main() -> std::io::Result<()> {
/// let path = Path::new("tests/data/test.fastq");
/// let mut reader = FastqReader::from_path(OptionPair::Single(path), 0, 0)?;
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
pub struct FastqReader<R: Read + Send> {
    inner: OptionPair<QReader<R>>,
    file_index: usize,
    reads_index: usize,
    // 批量读取
    batch_size: usize,

    // 新增：用于 Single 模式的状态
    single_mode: SingleReadMode,
    // 保留：用于“预读”的状态转换缓冲区
    read_ahead_buffer: Option<(Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>)>,
}

impl<R> FastqReader<R>
where
    R: Read + Send,
{
    /// Creates a new FastqReader with default capacity and batch size.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqkmer::{FastqReader, OptionPair};
    /// use std::fs::File;
    ///
    /// # fn main() -> std::io::Result<()> {
    /// let file = File::open("tests/data/test.fastq")?;
    /// let reader = FastqReader::new(OptionPair::Single(file), 0, 0);
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(readers: OptionPair<R>, file_index: usize, quality_score: i32) -> Self {
        Self::with_capacity(readers, file_index, BUFSIZE, quality_score, 30)
    }

    /// Creates a new FastqReader with specified capacity and batch size.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqkmer::{FastqReader, OptionPair};
    /// use std::fs::File;
    ///
    /// # fn main() -> std::io::Result<()> {
    /// let file = File::open("tests/data/test.fastq")?;
    /// let reader = FastqReader::with_capacity(OptionPair::Single(file), 0, 4096, 0, 50);
    /// # Ok(())
    /// # }
    /// ```
    pub fn with_capacity(
        readers: OptionPair<R>,
        file_index: usize,
        capacity: usize,
        quality_score: i32,
        batch_size: usize,
    ) -> Self {
        assert!(capacity >= 3);
        let inner = match readers {
            OptionPair::Single(reader) => {
                OptionPair::Single(QReader::with_capacity(reader, capacity, quality_score))
            }
            OptionPair::Pair(reader1, reader2) => OptionPair::Pair(
                QReader::with_capacity(reader1, capacity, quality_score),
                QReader::with_capacity(reader2, capacity, quality_score),
            ),
        };
        Self {
            inner,
            file_index,
            reads_index: 0,
            batch_size,
            single_mode: SingleReadMode::Unknown, // <-- 新增
            read_ahead_buffer: None,             // <-- 新增 (或保留)
        }
    }


    /// 辅助函数：从 header 字节中提取原始 ID 和基础 ID
    fn read_number_from_token(token: &str) -> Option<u8> {
        if token.is_empty() {
            return None;
        }

        if token.ends_with("/1") {
            return Some(1);
        }
        if token.ends_with("/2") {
            return Some(2);
        }

        let rest = token.as_bytes().get(1..);
        match token.as_bytes().first() {
            Some(b'1') if rest.map_or(true, |r| matches!(r.first(), Some(b':') | Some(b'/') | Some(b'#'))) => {
                Some(1)
            }
            Some(b'2') if rest.map_or(true, |r| matches!(r.first(), Some(b':') | Some(b'/') | Some(b'#'))) => {
                Some(2)
            }
            _ => None,
        }
    }

    fn get_ids_from_header(header_vec: &[u8]) -> (String, String, Option<u8>) {
        // 确保它能正确处理换行符
        if header_vec.is_empty() || header_vec[0] != b'@' {
            return (String::new(), String::new(), None);
        }

        let mut end = header_vec.len();
        if end > 0 && header_vec[end - 1] == b'\n' { end -= 1; }
        if end > 0 && header_vec[end - 1] == b'\r' { end -= 1; }
        
        if end <= 1 { return (String::new(), String::new(), None); }

        let s = unsafe {
            std::str::from_utf8_unchecked(&header_vec[1..end])
        };
        
        let first_space_index = s
            .find(|c: char| c.is_whitespace() || c == '\u{1}')
            .unwrap_or(s.len());
        
        let raw_id_str = &s[..first_space_index];
        let raw_id = raw_id_str.to_string();
        let base_id = trim_pair_info(raw_id_str);

        let mut read_number = Self::read_number_from_token(raw_id_str);
        if read_number.is_none() && first_space_index < s.len() {
            let rest = s[first_space_index..].trim_start();
            for token in rest.split_whitespace() {
                if let Some(num) = Self::read_number_from_token(token) {
                    read_number = Some(num);
                    break;
                }
            }
        }
        
        (raw_id, base_id, read_number)
    }


    fn create_seq_header(reader: &QReader<R>, file_index: usize, reads_index: usize) -> SeqHeader {
        let seq_id = unsafe {
            let s = std::str::from_utf8_unchecked(&reader.header[1..]);
            let first_space_index = s
                .find(|c: char| c.is_whitespace() || c == '\u{1}')
                .unwrap_or(s.len());

            // 直接从原始切片创建第一个单词的切片
            &s[..first_space_index]
        };
        SeqHeader {
            file_index,
            reads_index,
            format: SeqFormat::Fastq,
            id: trim_pair_info(seq_id),
        }
    }

    /// # Examples
    ///
    /// ```
    /// // 假设这些 struct 都在 crate 根部被 pub use
    /// use seqkmer::{FastqReader, OptionPair, Base, SeqHeader, SeqFormat};
    /// use std::io::Result;
    ///
    /// # fn main() -> Result<()> {
    ///
    /// // --- Test 1: Standard Paired-End (Two Files) ---
    /// let r1_data: &[u8] = b"@read1/1\nAAC\n+\nFFF\n";
    /// let r2_data: &[u8] = b"@read1/2\nTTG\n+\nIII\n";
    /// let readers_pair = OptionPair::Pair(r1_data, r2_data);
    /// let mut reader_pair = FastqReader::new(readers_pair, 0, 0);
    ///
    /// // 读取唯一的配对
    /// if let Some(seq_pair) = reader_pair.read_next()? {
    ///     assert_eq!(seq_pair.header.id, "read1"); // ID 被 trim_pair_info 修剪
    ///     match seq_pair.body {
    ///         OptionPair::Pair(r1, r2) => {
    ///             assert_eq!(r1, b"AAC");
    ///             assert_eq!(r2, b"TTG");
    ///         }
    ///         _ => panic!("Test 1 Failed: Expected Pair"),
    ///     }
    /// } else {
    ///     panic!("Test 1 Failed: Did not read any sequence");
    /// }
    /// // 检查文件是否结束
    /// assert!(reader_pair.read_next()?.is_none());
    ///
    ///
    /// // --- Test 2: Interleaved Paired-End (One File) ---
    /// let interleaved_data: &[u8] = b"@pairA/1\nGG\n+\n!!\n@pairA/2\nCC\n+\n!!\n@pairB/1\nT\n+\n#\n@pairB/2\nA\n+\n#\n";
    /// let mut reader_interleaved = FastqReader::new(OptionPair::Single(interleaved_data), 1, 0);
    ///
    /// // 读取第一个配对 (测试 Unknown -> Interleaved 状态转换)
    /// if let Some(seq_a) = reader_interleaved.read_next()? {
    ///     assert_eq!(seq_a.header.id, "pairA");
    ///     match seq_a.body {
    ///         OptionPair::Pair(r1, r2) => {
    ///             assert_eq!(r1, b"GG");
    ///             assert_eq!(r2, b"CC");
    ///         }
    ///         _ => panic!("Test 2 Failed: Expected Pair for seq_a"),
    ///     }
    /// } else {
    ///     panic!("Test 2 Failed: Did not read seq_a");
    /// }
    ///
    /// // 读取第二个配对 (测试 Interleaved 快速路径)
    /// if let Some(seq_b) = reader_interleaved.read_next()? {
    ///     assert_eq!(seq_b.header.id, "pairB");
    ///     match seq_b.body {
    ///         OptionPair::Pair(r1, r2) => {
    ///             assert_eq!(r1, b"T");
    ///             assert_eq!(r2, b"A");
    ///         }
    ///         _ => panic!("Test 2 Failed: Expected Pair for seq_b"),
    ///     }
    /// } else {
    ///     panic!("Test 2 Failed: Did not read seq_b");
    /// }
    ///
    /// // 检查文件是否结束
    /// assert!(reader_interleaved.read_next()?.is_none());
    ///
    ///
    /// // --- Test 3: Standard Single-End (One File) ---
    /// // ID 没有 /1 或 /2
    /// let single_data: &[u8] = b"@seq1\nATAT\n+\nFFFF\n@seq2_no_pair\nCGCG\n+\nIIII\n";
    /// let mut reader_single = FastqReader::new(OptionPair::Single(single_data), 2, 0);
    ///
    /// // 读取第一个 read (测试 Unknown -> SingleEnd 状态转换)
    /// if let Some(seq1) = reader_single.read_next()? {
    ///     assert_eq!(seq1.header.id, "seq1");
    ///     match seq1.body {
    ///         OptionPair::Single(r) => assert_eq!(r, b"ATAT"),
    ///         _ => panic!("Test 3 Failed: Expected Single for seq1"),
    ///     }
    /// } else {
    ///     panic!("Test 3 Failed: Did not read seq1");
    /// }
    ///
    /// // 读取第二个 read (测试 read_ahead_buffer 和 SingleEnd 快速路径)
    /// if let Some(seq2) = reader_single.read_next()? {
    ///     assert_eq!(seq2.header.id, "seq2_no_pair");
    ///     match seq2.body {
    ///         OptionPair::Single(r) => assert_eq!(r, b"CGCG"),
    ///         _ => panic!("Test 3 Failed: Expected Single for seq2"),
    ///     }
    /// } else {
    ///     panic!("Test 3 Failed: Did not read seq2");
    /// }
    ///
    /// // 检查文件是否结束
    /// assert!(reader_single.read_next()?.is_none());
    ///
    /// # Ok(())
    /// # }
    /// ```
    pub fn read_next(&mut self) -> Result<Option<Base<Vec<u8>>>> {
        match &mut self.inner {
            OptionPair::Single(reader) => {
                
                // --- 状态机 ---
                match self.single_mode {
                    
                    // --- 快速路径 1: Interleaved ---
                    SingleReadMode::Interleaved => {
                        // 我们已“锁定”为 Interleaved 模式。
                        // 假设文件格式是正确的，我们只管成对读取。
                        
                        // 读取 R1
                        if reader.read_next()?.is_none() {
                            return Ok(None); // 正常的文件结尾
                        }
                        self.reads_index += 1;
                        let seq_header = Self::create_seq_header(&reader, self.file_index, self.reads_index);
                        let r1_seq = reader.seq.to_owned();
                        
                        // 读取 R2
                        if reader.read_next()?.is_none() {
                            // 错误：文件不完整，R1 缺少 R2
                            return Err(std::io::Error::new(
                                std::io::ErrorKind::UnexpectedEof,
                                "Truncated interleaved FASTQ: R1 found but R2 is missing.",
                            ));
                        }
                        let r2_seq = reader.seq.to_owned();

                        Ok(Some(Base::new(
                            seq_header, // 使用 R1 的 header
                            OptionPair::Pair(r1_seq, r2_seq),
                        )))
                    }
                    
                    // --- 快速路径 2: SingleEnd ---
                    SingleReadMode::SingleEnd => {
                        // 我们已“锁定”为 SingleEnd 模式。
                        
                        // 首先检查缓冲区（仅在第一次转换时使用）
                        if let Some(buffered) = self.read_ahead_buffer.take() {
                             self.reads_index += 1;
                             let (_, base_id, _) = Self::get_ids_from_header(&buffered.0); // (buffered.0 is header)
                             let seq_header = SeqHeader {
                                file_index: self.file_index,
                                reads_index: self.reads_index,
                                format: SeqFormat::Fastq,
                                id: base_id,
                             };
                             return Ok(Some(Base::new(
                                seq_header,
                                OptionPair::Single(buffered.1), // (buffered.1 is seq)
                             )));
                        }

                        // 缓冲区为空，正常从文件读取
                        if reader.read_next()?.is_none() {
                            return Ok(None);
                        }
                        self.reads_index += 1;
                        let seq_header = Self::create_seq_header(&reader, self.file_index, self.reads_index);
                        Ok(Some(Base::new(
                            seq_header,
                            OptionPair::Single(reader.seq.to_owned()),
                        )))
                    }
                    
                    // --- 状态 3: Unknown (仅运行一次) ---
                    SingleReadMode::Unknown => {
                        // 这是第一次调用，执行昂贵的检测
                        
                        // 1. 读取 R1
                        if reader.read_next()?.is_none() {
                            return Ok(None); // 空文件
                        }
                        let r1_header_vec = reader.header.to_owned();
                        let r1_seq_vec = reader.seq.to_owned();
                        
                        // 2. 解析 R1
                        let (r1_raw_id, r1_base_id, r1_read_num) =
                            Self::get_ids_from_header(&r1_header_vec);
                        self.reads_index += 1;
                        let seq_header = SeqHeader {
                            file_index: self.file_index,
                            reads_index: self.reads_index,
                            format: SeqFormat::Fastq,
                            id: r1_base_id.clone(),
                        };

                        // 3. 尝试读取 R2
                        if reader.read_next()?.is_none() {
                            // 文件只有 1 条 read。它必须是 SingleEnd。
                            self.single_mode = SingleReadMode::SingleEnd;
                            return Ok(Some(Base::new(
                                seq_header,
                                OptionPair::Single(r1_seq_vec),
                            )));
                        }

                        // 4. 解析 R2 并检查
                        let (r2_raw_id, r2_base_id, r2_read_num) =
                            Self::get_ids_from_header(&reader.header);
                        let is_pair = r1_base_id == r2_base_id
                            && match (r1_read_num, r2_read_num) {
                                (Some(1), Some(2)) => true,
                                _ => r1_raw_id.ends_with("/1") && r2_raw_id.ends_with("/2"),
                            };

                        if is_pair {
                            // --- 检测为 Interleaved ---
                            self.single_mode = SingleReadMode::Interleaved;
                            Ok(Some(Base::new(
                                seq_header,
                                OptionPair::Pair(r1_seq_vec, reader.seq.to_owned()),
                            )))
                        } else {
                            // --- 检测为 SingleEnd ---
                            self.single_mode = SingleReadMode::SingleEnd;
                            
                            // 必须将 R2 存入缓冲区！
                            self.read_ahead_buffer = Some((
                                reader.header.to_owned(),
                                reader.seq.to_owned(),
                                reader.plus.to_owned(),
                                reader.quals.to_owned()
                            ));
                            
                            // 仅返回 R1
                            Ok(Some(Base::new(
                                seq_header,
                                OptionPair::Single(r1_seq_vec),
                            )))
                        }
                    }
                }
            }
            OptionPair::Pair(reader1, reader2) => {
                // --- 双文件模式逻辑不变 ---
                if reader1.read_next()?.is_none() {
                    return Ok(None);
                }
                if reader2.read_next()?.is_none() {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::UnexpectedEof,
                        "R2 file is shorter than R1 file.",
                    ));
                }

                self.reads_index += 1;
                let seq_header =
                    Self::create_seq_header(&reader1, self.file_index, self.reads_index);
                
                Ok(Some(Base::new(
                    seq_header,
                    OptionPair::Pair(reader1.seq.to_owned(), reader2.seq.to_owned()),
                )))
            }
        }
    }
}

impl FastqReader<Box<dyn Read + Send>> {
    /// Creates a new FastqReader from file paths.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqkmer::{FastqReader, OptionPair};
    /// use std::path::Path;
    ///
    /// # fn main() -> std::io::Result<()> {
    /// let path = Path::new("tests/data/test.fastq");
    /// let reader = FastqReader::from_path(OptionPair::Single(path), 0, 0)?;
    /// # Ok(())
    /// # }
    /// ```
    #[inline]
    pub fn from_path<P: AsRef<Path>>(
        paths: OptionPair<P>,
        file_index: usize,
        quality_score: i32,
    ) -> Result<Self> {
        let readers = paths.map(|path| dyn_reader(path))?;
        Ok(Self::new(readers, file_index, quality_score))
    }
}

impl<R> Reader for FastqReader<R>
where
    R: Read + Send,
{
    fn next(&mut self) -> Result<Option<Vec<Base<Vec<u8>>>>> {
        let seqs: Vec<Base<Vec<u8>>> = (0..self.batch_size)
            .filter_map(|_| self.read_next().transpose())
            .collect::<Result<Vec<_>>>()?;

        Ok(Some(seqs).filter(|v| !v.is_empty()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detects_interleaved_with_space_delimited_read_numbers() -> std::io::Result<()> {
        let data: &[u8] = b"@pairA 1:N:0:1\nAC\n+\n!!\n@pairA 2:N:0:1\nTG\n+\n!!\n";
        let mut reader = FastqReader::new(OptionPair::Single(data), 0, 0);

        let first_pair = reader.read_next()?;
        let base = first_pair.expect("expected a paired read");
        match base.body {
            OptionPair::Pair(r1, r2) => {
                assert_eq!(r1, b"AC");
                assert_eq!(r2, b"TG");
            }
            other => panic!("expected paired reads, got {:?}", other),
        }

        Ok(())
    }

    #[test]
    fn detects_interleaved_with_slash_suffixes() -> std::io::Result<()> {
        let data: &[u8] = b"@pair1/1\nAC\n+\n!!\n@pair1/2\nTG\n+\n!!\n";
        let mut reader = FastqReader::new(OptionPair::Single(data), 0, 0);

        let base = reader
            .read_next()?
            .expect("expected a paired read in interleaved mode");
        assert_eq!(base.header.id, "pair1");
        match base.body {
            OptionPair::Pair(r1, r2) => {
                assert_eq!(r1, b"AC");
                assert_eq!(r2, b"TG");
            }
            other => panic!("expected paired reads, got {:?}", other),
        }

        Ok(())
    }

    #[test]
    fn single_end_fastq_remains_single() -> std::io::Result<()> {
        let data: &[u8] =
            b"@seq1\nAC\n+\n!!\n@seq2\nTT\n+\n!!\n";
        let mut reader = FastqReader::new(OptionPair::Single(data), 0, 0);

        let first = reader
            .read_next()?
            .expect("expected first single-end read");
        assert_eq!(first.header.id, "seq1");
        match first.body {
            OptionPair::Single(seq) => assert_eq!(seq, b"AC"),
            other => panic!("expected single body, got {:?}", other),
        }

        let second = reader
            .read_next()?
            .expect("expected second single-end read");
        assert_eq!(second.header.id, "seq2");
        match second.body {
            OptionPair::Single(seq) => assert_eq!(seq, b"TT"),
            other => panic!("expected single body, got {:?}", other),
        }

        assert!(reader.read_next()?.is_none());

        Ok(())
    }

    #[test]
    fn interleaved_reader_errors_when_r2_missing() {
        let data: &[u8] = b"@pair/1\nAA\n+\n!!\n@pair/2\nTT\n+\n!!\n@broken/1\nGG\n+\n!!\n";
        let mut reader = FastqReader::new(OptionPair::Single(data), 0, 0);

        // First pair should succeed.
        assert!(reader.read_next().unwrap().is_some());

        // Second pair is truncated: expect an UnexpectedEof error.
        let err = reader.read_next().unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::UnexpectedEof);
    }

    #[test]
    fn paired_files_mode_requires_balanced_lengths() {
        let r1: &[u8] = b"@r1/1\nAC\n+\n!!\n";
        let r2: &[u8] = b"";
        let mut reader = FastqReader::new(OptionPair::Pair(r1, r2), 0, 0);

        let err = reader.read_next().unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::UnexpectedEof);
    }

    #[test]
    fn quality_threshold_masks_low_quality_bases() -> std::io::Result<()> {
        let data: &[u8] = b"@seq\nACGT\n+\n!I!I\n";
        let mut reader = FastqReader::new(OptionPair::Single(data), 0, 10);

        let base = reader
            .read_next()?
            .expect("expected single read after quality filtering");

        match base.body {
            OptionPair::Single(seq) => assert_eq!(seq, b"xCxT"),
            other => panic!("expected single body, got {:?}", other),
        }

        Ok(())
    }
}
