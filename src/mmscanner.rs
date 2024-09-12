// kraken 2 使用的是murmur_hash3 算法的 fmix64作为 hash
use crate::feat::Meros;
use crate::feat::{canonical_representation, char_to_value, fmix64 as murmur_hash3};
use crate::Base;
use crate::OptionPair;
use crate::BITS_PER_CHAR;
use std::collections::VecDeque;

/// Converts an l-mer to a candidate l-mer using the given Meros configuration.
#[inline]
fn to_candidate_lmer(meros: &Meros, lmer: u64) -> u64 {
    let mut canonical_lmer = canonical_representation(lmer, meros.l_mer);
    if meros.spaced_seed_mask > 0 {
        canonical_lmer &= meros.spaced_seed_mask;
    }
    canonical_lmer ^ meros.toggle_mask
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_to_candidate_lmer() {
        let meros = Meros::new(11, 3, Some(0), None, None);
        let lmer = 0b11001100110011001100u64;
        let candidate = to_candidate_lmer(&meros, lmer);
        // println!("Candidate l-mer: {:b}", candidate);
        // 在这里添加断言来验证结果
        assert_eq!(candidate, 0b11110u64);
    }
}

/// Represents data for a minimizer.
///
/// # Examples
///
/// ```
/// use seqkmer::MinimizerData;
///
/// let data = MinimizerData::new(0b1100110011u64, 5);
/// println!("Position: {}, Candidate l-mer: {:b}", data.pos, data.candidate_lmer);
/// ```
#[derive(Debug)]
pub struct MinimizerData {
    pub pos: usize,
    pub candidate_lmer: u64,
}

impl MinimizerData {
    /// Creates a new MinimizerData instance.
    pub fn new(candidate_lmer: u64, pos: usize) -> Self {
        Self {
            candidate_lmer,
            pos,
        }
    }
}

/// Represents a window for finding minimizers.
///
/// # Examples
///
/// ```
/// use seqkmer::MinimizerWindow;
///
/// let mut window = MinimizerWindow::new(5);
/// let minimizer = window.next(0b1100110011u64);
/// println!("Minimizer: {:?}", minimizer);
/// ```
pub struct MinimizerWindow {
    queue: VecDeque<MinimizerData>,
    queue_pos: usize,
    /// 窗口队列的大小
    capacity: usize,
    /// 队列计数
    count: usize,
}

impl MinimizerWindow {
    /// Creates a new MinimizerWindow with the given capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            queue: VecDeque::with_capacity(capacity),
            capacity,
            count: 0,
            queue_pos: 0,
        }
    }

    /// Processes the next candidate l-mer and returns the minimizer if it has changed.
    #[inline]
    pub fn next(&mut self, candidate_lmer: u64) -> Option<u64> {
        // 无需比较，直接返回
        if self.capacity == 1 {
            return Some(candidate_lmer);
        }

        let data = MinimizerData::new(candidate_lmer, self.count);

        // 移除队列中所有比当前元素大的元素的索引
        // 因为它们不可能是当前窗口的最小值
        while let Some(m_data) = self.queue.back() {
            if m_data.candidate_lmer > candidate_lmer {
                self.queue.pop_back();
            } else {
                break;
            }
        }
        let mut changed = false;

        if (self.queue.is_empty() && self.count >= self.capacity) || self.count == self.capacity {
            changed = true
        }
        // 将当前元素的索引添加到队列
        self.queue.push_back(data);

        while !self.queue.is_empty()
            && self.queue.front().map_or(false, |front| {
                self.count >= self.capacity && front.pos < self.count - self.capacity
            })
        {
            self.queue.pop_front();
            changed = true;
        }

        self.count += 1;
        if changed {
            self.queue.front().map(|front| front.candidate_lmer)
        } else {
            None
        }
    }

    /// Clears the window.
    fn clear(&mut self) {
        self.count = 0;
        self.queue_pos = 0;
        self.queue.clear();
    }
}

/// Represents a cursor for processing l-mers.
///
/// # Examples
///
/// ```
/// use seqkmer::{Cursor, Meros};
///
/// let meros = Meros::new(11, 3, Some(0), None, None);
/// let cursor = Cursor::new(&meros);
/// ```
#[derive(Clone, Copy)]
pub struct Cursor {
    pos: usize,
    capacity: usize,
    value: u64,
    mask: u64,
}

impl Cursor {
    /// Creates a new Cursor with the given capacity and mask.
    pub fn new(meros: &Meros) -> Self {
        Self {
            pos: 0,
            value: 0,
            capacity: meros.l_mer,
            mask: meros.mask,
        }
    }

    /// Processes the next item and returns the l-mer if the capacity is reached.
    fn next_lmer(&mut self, item: u64) -> Option<u64> {
        self.value = ((self.value << BITS_PER_CHAR) | item) & self.mask;
        // 更新当前位置
        self.pos += 1;
        // 检查是否达到了容量
        if self.pos >= self.capacity {
            return Some(self.value);
        }
        None
    }

    /// Clears the cursor.
    #[inline]
    fn clear(&mut self) {
        self.pos = 0;
        self.value = 0;
    }
}

/// Iterator for finding minimizers in a sequence.
///
/// # Examples
///
/// ```
/// use seqkmer::{MinimizerIterator, Meros, Cursor, MinimizerWindow};
///
/// let meros = Meros::new(11, 3, Some(0), None, None);
/// let cursor = Cursor::new(&meros);
/// let window = MinimizerWindow::new(meros.window_size());
/// let mut iter = MinimizerIterator::new(b"ATCGATCGATCG", cursor, window, &meros);
///
/// for (pos, minimizer) in iter {
///     println!("Position: {}, Minimizer: {:b}", pos, minimizer);
/// }
/// ```
pub struct MinimizerIterator<'a> {
    cursor: Cursor,
    window: MinimizerWindow,
    seq: &'a [u8],
    meros: &'a Meros,
    pos: usize,
    end: usize,
    pub size: usize,
}

impl<'a> MinimizerIterator<'a> {
    /// Creates a new MinimizerIterator.
    pub fn new(seq: &'a [u8], cursor: Cursor, window: MinimizerWindow, meros: &'a Meros) -> Self {
        MinimizerIterator {
            cursor,
            window,
            seq,
            meros,
            pos: 0,
            size: 0,
            end: seq.len(),
        }
    }

    /// Clears the internal state of the iterator.
    fn clear_state(&mut self) {
        self.cursor.clear();
        self.window.clear();
    }

    /// Returns the size of the sequence being processed.
    pub fn seq_size(&self) -> usize {
        self.end
    }
}

impl<'a> Iterator for MinimizerIterator<'a> {
    type Item = (usize, u64);

    fn next(&mut self) -> Option<Self::Item> {
        while self.pos < self.end {
            let ch = self.seq[self.pos];
            self.pos += 1;
            if ch == b'\n' || ch == b'\r' {
                continue;
            } else {
                let data = match char_to_value(ch) {
                    Some(code) => self.cursor.next_lmer(code).and_then(|lmer| {
                        let candidate_lmer = to_candidate_lmer(&self.meros, lmer);
                        self.window
                            .next(candidate_lmer)
                            .map(|minimizer| murmur_hash3(minimizer ^ self.meros.toggle_mask))
                    }),
                    None => {
                        self.clear_state();
                        None
                    }
                };
                if data.is_some() {
                    self.size += 1;
                    return Some((self.size, data.unwrap()));
                }
            }
        }
        None
    }
}

impl<'a> Base<MinimizerIterator<'a>> {
    /// Returns the size of the sequence as a string.
    pub fn seq_size_str(&self) -> OptionPair<String> {
        self.body.apply(|m_iter| m_iter.seq_size().to_string())
    }

    /// Formats the size of the sequence as a string.
    pub fn fmt_seq_size(&self) -> String {
        self.body
            .reduce_str("|", |m_iter| m_iter.seq_size().to_string())
    }

    /// Formats the size of the minimizer iterator as a string.
    pub fn fmt_size(&self) -> String {
        self.body.reduce_str("|", |m_iter| m_iter.size.to_string())
    }

    /// Folds the minimizer iterator into a vector of a specified type.
    pub fn fold<F, T>(&mut self, mut f: F) -> Vec<T>
    where
        F: FnMut(&mut Vec<T>, &mut MinimizerIterator<'a>, usize) -> usize,
        T: Clone,
    {
        let mut init = Vec::new();
        match &mut self.body {
            OptionPair::Single(m_iter) => {
                f(&mut init, m_iter, 0);
            }
            OptionPair::Pair(m_iter1, m_iter2) => {
                let offset = f(&mut init, m_iter1, 0);
                f(&mut init, m_iter2, offset);
            }
        }
        init
    }

    /// Returns the range of the minimizer iterator.
    pub fn range(&self) -> OptionPair<(usize, usize)> {
        match &self.body {
            OptionPair::Single(m_iter) => OptionPair::Single((0, m_iter.size)),
            OptionPair::Pair(m_iter1, m_iter2) => {
                let size1 = m_iter1.size;
                OptionPair::Pair((0, size1), (size1, m_iter2.size + size1))
            }
        }
    }
}

/// Scans a sequence and returns a MinimizerIterator.
///
/// # Examples
///
/// ```
/// use seqkmer::{scan_sequence, Base, Meros, OptionPair, SeqHeader, SeqFormat};
///
/// let header = SeqHeader {
///     id: "test".to_string(),
///     file_index: 0,
///     reads_index: 0,
///     format: SeqFormat::Fasta,
/// };
/// let seq = Base::new(header, OptionPair::Single(b"ATCGATCGATCG".to_vec()));
/// let meros = Meros::new(11, 3, Some(0), None, None);
/// let mut minimizer_iter = scan_sequence(&seq, &meros);
///
/// if let OptionPair::Single(ref mut iter) = minimizer_iter.body {
///     for minimizer in iter {
///         println!("Minimizer: {:?}", minimizer);
///     }
/// }
/// ```
pub fn scan_sequence<'a>(
    sequence: &'a Base<Vec<u8>>,
    meros: &'a Meros,
) -> Base<MinimizerIterator<'a>> {
    let func = |seq: &'a Vec<u8>| {
        let cursor = Cursor::new(meros);
        let window = MinimizerWindow::new(meros.window_size());
        MinimizerIterator::new(seq, cursor, window, meros)
    };

    match &sequence.body {
        OptionPair::Pair(seq1, seq2) => Base::new(
            sequence.header.clone(),
            OptionPair::Pair(func(&seq1), func(&seq2)),
        ),
        OptionPair::Single(seq1) => {
            Base::new(sequence.header.clone(), OptionPair::Single(func(&seq1)))
        }
    }
}
