# Seqkmer

Seqkmer is a Rust-based tool for processing and analyzing DNA sequences using k-mers.

## Description

Seqkmer is designed to efficiently handle FASTA files and perform k-mer analysis on DNA sequences. It provides functionality for reading FASTA files, generating k-mers, and potentially other sequence analysis tasks.

## Features

- FASTA file parsing
- K-mer generation


## Usage

```rust
use seqkmer::{FastxReader, Reader};
use std::path::Path;

fn main() -> std::io::Result<()> {
    // Create a FastxReader for a single FASTA file
    let path = Path::new("path/to/your/fasta_file.fa");
    let mut reader = FastxReader::from_paths(
        seqkmer::OptionPair::Single(path),
        0, // file index
        0  // quality score (not used for FASTA)
    )?;

    // Read sequences
    while let Some(sequences) = reader.next()? {
        for sequence in sequences {
            println!("Sequence ID: {}", sequence.header.id);
            println!("Sequence length: {}", sequence.body.single().unwrap().len());
            // Process the sequence as needed
        }
    }

    Ok(())
}
```
