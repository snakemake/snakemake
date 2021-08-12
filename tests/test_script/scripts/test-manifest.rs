//! This is a regular crate doc comment, but it also contains a partial
//! Cargo manifest.  Note the use of a *fenced* code block, and the
//! `cargo` "language".
//!
//! ```cargo
//! [dependencies]
//! csv = "1.1"
//! serde = { version = "1.0", features = ["derive"] }
//! ```


use std::error::Error;
use std::io::{BufWriter, Write};
use std::fs::File;
use serde::Deserialize;


static BED: &[u8] = b"chrom1	1	15	foo	454	-
chrom1	40	45	bar	2	+
chrom2	4	45	baz	2	-
";

#[derive(Debug, Deserialize)]
struct BedRecord {
    chrom: String,
    start: u64,
    end: u64,
    name: Option<String>,
    score: Option<u16>,
    strand: Option<char>,
}

fn main() -> Result<(), Box<dyn Error>> {
    snakemake.redirect_stderr(&snakemake.log[0])?;
    let f_out = File::create(&snakemake.output[0])?;

    let mut ostream = BufWriter::new(f_out);
    println!("Loaded");

    let keep_strand = match &snakemake.params.keep {
        s if s.len() == 1 => Some(s.chars().next().unwrap() as char),
        _ => None,
    };

    println!("Reading BED file...");
    let mut rdr = csv::ReaderBuilder::new().has_headers(false).delimiter(b'\t').from_reader(BED);
    for result in rdr.deserialize() {
        // Notice that we need to provide a type hint for automatic
        // deserialization.
        let record: BedRecord = result?;
        let l = record.end - record.start;
        if record.strand == keep_strand {
            write!(&mut ostream, "{}\t{}\n", record.chrom, l)?;
        }
    }
    println!("Output written to {}", &snakemake.output[0]);
    Ok(())
}
