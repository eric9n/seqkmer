# Seqkmer

Seqkmer is a Rust library for high-throughput sequence IO and k-mer based analyses. It provides fast readers for FASTA/FASTQ (including gzipped streams), k-mer minimizer scanning, and utilities to parallelise bulk sequence processing.

## Highlights

- **Universal FASTX readers**: Seamlessly handle FASTA, FASTQ, interleaved paired-end, and dual-file paired-end datasets through a unified API. Automatic format detection and transparent gzip support are included.
- **Quality-aware FASTQ parsing**: Optional quality-score thresholds to soft-mask low-quality bases while preserving original sequence layout.
- **Buffered & streaming modes**: Choose between streaming (`FastaReader`, `FastqReader`) or buffered variants (`BufferFastaReader`) depending on your throughput/memory trade-offs.
- **Minimizer-based k-mer scanning**: The `mmscanner` module exposes `scan_sequence` and `MinimizerIterator` for fast k-mer/minimizer enumeration with configurable windows.
- **Parallel orchestration**: Utilities in `parallel` coordinate multi-threaded reading and processing pipelines using scoped thread pools.

## Getting Started

Add Seqkmer to your project:

```bash
cargo add seqkmer
```

### Reading FASTA or FASTQ

```rust
use seqkmer::{FastxReader, OptionPair, Reader};
use std::path::Path;

fn main() -> std::io::Result<()> {
    // Single FASTQ file (auto-detects FASTA vs FASTQ and gzip)
    let path = Path::new("tests/data/test.fastq");
    let mut reader = FastxReader::from_paths(OptionPair::Single(path), 0, 18)?;

    while let Some(batch) = reader.next()? {
        for entry in batch {
            println!(
                "[{}] {} (len={})",
                entry.header.format as u8,
                entry.header.id,
                entry.body.single().unwrap().len()
            );
        }
    }

    Ok(())
}
```

For paired-end data, provide a pair of paths. Interleaved FASTQ is detected automatically; separate R1/R2 files are also supported:

```rust
let paths = OptionPair::Pair(
    Path::new("reads_R1.fastq"),
    Path::new("reads_R2.fastq"),
);
let mut reader = FastxReader::from_paths(paths, 0, 0)?;
```

### K-mer Minimizer Scanning

```rust
use seqkmer::{scan_sequence, Meros, MinimizerIterator};
use seqkmer::reader::Reader;

fn main() -> std::io::Result<()> {
    let meros = Meros::new(15, 5, Some(0), None, None); // (k, window, seed, min, max)
    let mut reader = seqkmer::FastaReader::from_path("tests/data/test.fasta", 0)?;

    while let Some(batch) = reader.next()? {
        for base in batch {
            let mut minimizers: Vec<_> = scan_sequence(&base, &meros).collect();
            println!("{} -> {} minimizers", base.header.id, minimizers.len());
        }
    }

    Ok(())
}
```

### Parallel Pipelines

Use `read_parallel` when you need to map a function across batches using multiple threads:

```rust
use seqkmer::{read_parallel, FastaReader, Meros, ParallelResult, Reader};

fn main() -> std::io::Result<()> {
    let meros = Meros::new(11, 3, Some(0), None, None);
    let mut reader = FastaReader::from_path("tests/data/test.fasta", 0)?;

    read_parallel(
        &mut reader,
        4, // threads
        &meros,
        |seqs| seqs.len(), // work: count sequences per batch
        |result: &mut ParallelResult<usize>| {
            let mut total = 0;
            while let Some(count) = result.next() {
                total += count.unwrap();
            }
            println!("processed {} batches", total);
        },
    )?;

    Ok(())
}
```

## Feature Overview

| Module          | Purpose                                                                 |
|-----------------|-------------------------------------------------------------------------|
| `fasta`         | FASTA readers (streaming + buffered)                                     |
| `fastq`         | FASTQ reader with automatic interleaved detection and quality masking    |
| `fastx`         | Format-agnostic wrapper over FASTA/FASTQ readers                         |
| `reader`        | Misc IO utilities (gzip detection, trim helpers, file format detection) |
| `parallel`      | Threaded reader orchestration using scoped thread pools                  |
| `mmscanner`     | Minimizer scanning over DNA sequences                                    |
| `feat`          | K-mer feature helper types (`Meros`, constants)                          |
| `utils::OptionPair` | Helper enum for representing single vs paired resources           |

## Testing

All functionality is covered by unit and doc tests. Run the full suite with:

```bash
cargo test
```

## License

Seqkmer is distributed under the terms of the MIT License.
