use crate::feat::Meros;
use crate::mmscanner::scan_sequence;
use crate::reader::detect_file_format;
use crate::reader::Reader;
use crate::seq::{Base, SeqFormat};
use crate::MinimizerIterator;
use crate::{FastaReader, FastqReader};
use crossbeam_channel::{bounded, Receiver};
use scoped_threadpool::Pool;
use std::collections::HashMap;
use std::io::Result;
use std::sync::Arc;

/// A wrapper for parallel processing items.
///
/// # Examples
///
/// ```
/// use seqkmer::ParallelItem;
///
/// let item = ParallelItem(42);
/// assert_eq!(item.unwrap(), 42);
/// ```
pub struct ParallelItem<P>(pub P);

impl<P> ParallelItem<P> {
    /// Unwraps the inner value.
    pub fn unwrap(self) -> P {
        self.0
    }
}

/// Represents the result of a parallel operation.
pub struct ParallelResult<P>
where
    P: Send,
{
    recv: Receiver<P>,
}

impl<P> ParallelResult<P>
where
    P: Send,
{
    /// Retrieves the next item from the parallel result.
    #[inline]
    pub fn next(&mut self) -> Option<ParallelItem<P>> {
        self.recv.recv().ok().map(ParallelItem)
    }
}

/// Creates a reader based on the file format.
///
/// # Examples
///
/// ```
/// use seqkmer::create_reader;
/// use std::path::Path;
///
/// # fn main() -> std::io::Result<()> {
/// let file_path = Path::new("tests/data/test.fasta").to_str().unwrap().to_string();
/// let reader = create_reader(&[file_path], 0, 0)?;
/// # Ok(())
/// # }
/// ```
pub fn create_reader(
    file_pair: &[String],
    file_index: usize,
    score: i32,
) -> Result<Box<dyn Reader + Send>> {
    // let mut files_iter = file_pair.iter();
    let paths = crate::OptionPair::from_slice(file_pair);

    match detect_file_format(&file_pair[0])? {
        SeqFormat::Fastq => Ok(Box::new(FastqReader::from_path(paths, file_index, score)?)),
        SeqFormat::Fasta => Ok(Box::new(FastaReader::from_path(&file_pair[0], file_index)?)),
    }
}

/// Performs parallel reading and processing of sequences.
///
/// # Examples
///
/// ```
/// use seqkmer::{read_parallel, FastaReader, Meros, Base, MinimizerIterator};
/// use std::path::Path;
///
/// # fn main() -> std::io::Result<()> {
/// let path = Path::new("tests/data/test.fasta");
/// let mut reader = FastaReader::from_path(path, 0)?;
/// let meros = Meros::new(11, 3, Some(0), None, None);
///
/// let work = |seqs: &mut Vec<Base<MinimizerIterator>>| {
///     // Process sequences
///     seqs.len()
/// };
///
/// let func = |result: &mut seqkmer::ParallelResult<usize>| {
///     let mut total = 0;
///     while let Some(count) = result.next() {
///         total += count.unwrap();
///     }
///     total
/// };
///
/// let total = read_parallel(&mut reader, 4, &meros, work, func)?;
/// println!("Total sequences processed: {:?}", total);
/// # Ok(())
/// # }
/// ```
pub fn read_parallel<R, W, O, F, Out>(
    reader: &mut R,
    n_threads: usize,
    meros: &Meros,
    work: W,
    func: F,
) -> Result<()>
where
    R: Reader,
    O: Send,
    Out: Send + Default,
    W: Send + Sync + Fn(&mut Vec<Base<MinimizerIterator>>) -> O,
    F: FnOnce(&mut ParallelResult<O>) -> Out + Send,
{
    assert!(n_threads > 2);
    let buffer_len = n_threads + 2;
    let (sender, receiver) = bounded::<Vec<Base<Vec<u8>>>>(buffer_len);
    let (done_send, done_recv) = bounded::<O>(buffer_len);
    let receiver = Arc::new(receiver); // 使用 Arc 来共享 receiver
    let done_send = Arc::new(done_send);
    let mut pool = Pool::new(n_threads as u32);

    let mut parallel_result = ParallelResult { recv: done_recv };

    pool.scoped(|pool_scope| {
        // 生产者线程
        pool_scope.execute(move || {
            while let Ok(Some(seqs)) = reader.next() {
                sender.send(seqs).expect("Failed to send sequences");
            }
        });

        // 消费者线程
        for _ in 0..n_threads - 2 {
            let receiver = Arc::clone(&receiver);
            let work = &work;
            let done_send = Arc::clone(&done_send);
            pool_scope.execute(move || {
                while let Ok(mut seqs) = receiver.recv() {
                    let mut markers: Vec<Base<MinimizerIterator<'_>>> = seqs
                        .iter_mut()
                        .map(|seq| scan_sequence(seq, &meros))
                        .collect();
                    let output = work(&mut markers);
                    done_send.send(output).expect("Failed to send outputs");
                }
            });
        }

        // 引用计数减掉一个,这样都子线程结束时, done_send还能完全释放
        drop(done_send);
        pool_scope.execute(move || {
            let _ = func(&mut parallel_result);
        });

        pool_scope.join_all();
    });

    Ok(())
}

/// Performs parallel reading and processing of buffered data.
///
/// # Examples
///
/// ```
/// use seqkmer::{buffer_read_parallel, FastaReader};
/// use std::path::Path;
/// use std::fs::File;
///
/// # fn main() -> std::io::Result<()> {
/// let path = Path::new("tests/data/test.fasta");
/// let mut file = File::open(path)?;
///
/// let work = |data: Vec<u8>| {
///     // Process data
///     data.len()
/// };
///
/// let func = |result: &mut seqkmer::ParallelResult<usize>| {
///     let mut total = 0;
///     while let Some(count) = result.next() {
///         total += count.unwrap();
///     }
///     total
/// };
///
/// let _ = buffer_read_parallel(&mut file, 4, 1024, work, func)?;
/// # Ok(())
/// # }
/// ```
pub fn buffer_read_parallel<R, D, W, O, F, Out>(
    reader: &mut R,
    n_threads: usize,
    buffer_size: usize,
    work: W,
    func: F,
) -> Result<()>
where
    D: Send + Sized + Sync + Clone,
    R: std::io::Read + Send,
    O: Send,
    Out: Send + Default,
    W: Send + Sync + Fn(Vec<D>) -> O,
    F: FnOnce(&mut ParallelResult<O>) -> Out + Send,
{
    assert!(n_threads > 2);
    let buffer_len = n_threads + 2;
    let (sender, receiver) = bounded::<Vec<D>>(buffer_len);
    let (done_send, done_recv) = bounded::<O>(buffer_len);
    let receiver = Arc::new(receiver); // 使用 Arc 来共享 receiver
    let done_send = Arc::new(done_send);
    let mut pool = Pool::new(n_threads as u32);

    let slot_size = std::mem::size_of::<D>();
    let mut parallel_result = ParallelResult { recv: done_recv };

    pool.scoped(|pool_scope| {
        // 生产者线程
        pool_scope.execute(move || {
            let mut batch_buffer = vec![0u8; slot_size * buffer_size];

            while let Ok(bytes_read) = reader.read(&mut batch_buffer) {
                if bytes_read == 0 {
                    break;
                } // 文件末尾

                let slots_in_batch = bytes_read / slot_size;
                let slots = unsafe {
                    std::slice::from_raw_parts(batch_buffer.as_ptr() as *const D, slots_in_batch)
                };
                sender
                    .send(slots.to_vec())
                    .expect("Failed to send sequences");
            }
        });

        // 消费者线程
        for _ in 0..n_threads - 2 {
            let receiver = Arc::clone(&receiver);
            let work = &work;
            let done_send = Arc::clone(&done_send);
            pool_scope.execute(move || {
                while let Ok(seqs) = receiver.recv() {
                    let output = work(seqs);
                    done_send.send(output).expect("Failed to send outputs");
                }
            });
        }

        // 引用计数减掉一个,这样都子线程结束时, done_send还能完全释放
        drop(done_send);
        pool_scope.execute(move || {
            let _ = func(&mut parallel_result);
        });

        pool_scope.join_all();
    });

    Ok(())
}

/// Performs parallel processing on a HashMap.
///
/// # Examples
///
/// ```
/// use seqkmer::buffer_map_parallel;
/// use std::collections::HashMap;
///
/// # fn main() -> std::io::Result<()> {
/// let mut map = HashMap::new();
/// map.insert(1, vec![1, 2, 3]);
/// map.insert(2, vec![4, 5, 6]);
///
/// let work = |(key, value): (&u32, &Vec<i32>)| {
///     // Process data
///     value.iter().sum::<i32>()
/// };
///
/// let func = |result: &mut seqkmer::ParallelResult<i32>| {
///     let mut total = 0;
///     while let Some(sum) = result.next() {
///         total += sum.unwrap();
///     }
///     total
/// };
///
/// let total = buffer_map_parallel(&map, 4, work, func)?;
/// println!("Total sum: {:?}", total);
/// # Ok(())
/// # }
/// ```
pub fn buffer_map_parallel<D, W, O, F, Out>(
    map: &HashMap<u32, Vec<D>>,
    n_threads: usize,
    work: W,
    func: F,
) -> Result<()>
where
    D: Send + Sized + Sync,
    O: Send,
    Out: Send + Default,
    W: Send + Sync + Fn((&u32, &Vec<D>)) -> O,
    F: FnOnce(&mut ParallelResult<O>) -> Out + Send,
{
    assert!(n_threads > 2);
    let buffer_len = n_threads + 2;
    let (sender, receiver) = bounded::<(&u32, &Vec<D>)>(buffer_len);
    let (done_send, done_recv) = bounded::<O>(buffer_len);
    let receiver = Arc::new(receiver); // 使用 Arc 来共享 receiver
    let done_send = Arc::new(done_send);
    let mut pool = Pool::new(n_threads as u32);

    let mut parallel_result = ParallelResult { recv: done_recv };

    pool.scoped(|pool_scope| {
        // 生产者线程
        pool_scope.execute(move || {
            for entry in map {
                sender.send(entry).expect("Failed to send sequences");
            }
        });

        // 消费者线程
        for _ in 0..n_threads - 2 {
            let receiver = Arc::clone(&receiver);
            let work = &work;
            let done_send = Arc::clone(&done_send);
            pool_scope.execute(move || {
                while let Ok(seqs) = receiver.recv() {
                    let output = work(seqs);
                    done_send.send(output).expect("Failed to send outputs");
                }
            });
        }

        // 引用计数减掉一个,这样都子线程结束时, done_send还能完全释放
        drop(done_send);
        pool_scope.execute(move || {
            let _ = func(&mut parallel_result);
        });

        pool_scope.join_all();
    });

    Ok(())
}
