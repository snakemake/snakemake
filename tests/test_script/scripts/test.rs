use std::io::Write;
println!("Rust script executing!");

let snakemake = Snakemake::load()?;
dbg!(&snakemake);
assert_eq!(snakemake.config["test"], Value::Bool(true));
assert_eq!(snakemake.params["integer"], Value::I64(123));
assert_eq!(snakemake.output["0"], "rust.out");
assert_eq!(snakemake.input["0"], "test.in");
assert_eq!(snakemake.input["named_input"], "test.in");

let mut f = std::fs::File::create(&snakemake.output["0"])?;
write!(&mut f, "Rust test succeded!")?;