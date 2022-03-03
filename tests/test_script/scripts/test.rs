use std::io::Write;
println!("Rust script executing!");

assert_eq!(snakemake.config.test, true);
assert_eq!(snakemake.params.integer, 123);
assert_eq!(snakemake.output[0], "rust.out");
assert_eq!(snakemake.input[0], "test.in");
assert_eq!(snakemake.input.named_input, "test.in");
for (idx, val) in (&snakemake.input).into_iter().enumerate() {
    dbg!(idx, &val);
}

let input = &snakemake.input;
for value in input {
    dbg!(value);
}

let mut f = std::fs::File::create(&snakemake.output[0])?;
write!(&mut f, "Rust test succeded!")?;
