use clap::Parser;
use rayon::ThreadPoolBuilder;
use single_rust::io::FileScope;
use std::fs::read_to_string;
use std::path::Path;

pub(crate) mod hasher;

#[derive(Debug, Parser)]
struct Args {
    #[clap(short, long)]
    input: String,

    #[clap(long)]
    intersect: String,

    #[clap(short, long)]
    output: String,

    #[clap(short,long)]
    jobs: usize,

    #[clap(short,long)]
    seed: i64
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();
    let intersect = read_intersection_file(args.intersect)?;
    let adata = single_rust::io::read_h5ad(args.input, FileScope::Read, false)?;
    let thread_pool = ThreadPoolBuilder::new().num_threads(args.jobs).build()?;
    thread_pool.install(|| {
        hasher::hash_cells(adata, intersect, args.seed, &args.output)
    })
}

fn read_intersection_file<P: AsRef<Path>>(path: P) -> anyhow::Result<Vec<String>> {
    Ok(read_to_string(path)?.lines().map(String::from).collect())
}
