use anyhow::{anyhow, bail, Context, Result};
use clap::{Parser, Subcommand};
use cmd_lib::run_fun;
use semver::{Version, VersionReq};
use std::env;
use std::ffi::OsStr;
use std::fs;
use std::path::PathBuf;
use std::process::Command;
use which::which;

use serde::{Deserialize, Serialize};
use serde_json::json;

mod utils;
use utils::prog_utils::*;

#[derive(Debug, Subcommand)]
enum Commands {
    /// build the splici index
    #[clap(arg_required_else_help = true)]
    Build {
        /// reference genome
        #[clap(short, long, value_parser)]
        fasta: PathBuf,

        /// reference GTF file
        #[clap(short, long, value_parser)]
        gtf: PathBuf,

        /// the target read length the index will be built for
        #[clap(short, long, value_parser)]
        rlen: u32,

        /// path to output directory (will be created if it doesn't exist)
        #[clap(short, long, value_parser)]
        output: PathBuf,

        /// path to FASTA file with extra spliced sequence to add to the index
        #[clap(short, long, value_parser)]
        spliced: Option<PathBuf>,

        /// path to FASTA file with extra unspliced sequence to add to the index
        #[clap(short, long, value_parser)]
        unspliced: Option<PathBuf>,

        /// deduplicate identical sequences inside the R script when building the splici reference
        #[clap(short = 'd', long = "dedup", action)]
        dedup: bool,

        /// if this flag is passed, build the sparse rather than dense index for mapping
        #[clap(short = 'p', long = "sparse", action)]
        sparse: bool,

        /// number of threads to use when running [default: min(16, num cores)]"
        #[clap(short, long, default_value_t = 16, value_parser)]
        threads: u32,
    },
    /// quantify a sample
    #[clap(arg_required_else_help = true)]
    Quant {
        /// path to index
        #[clap(short, long, value_parser)]
        index: String,
    },
}

/// simplifying alevin-fry workflows
#[derive(Debug, Parser)]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}


fn main() -> anyhow::Result<()> {
    // gather information about the required 
    // programs.
    let rp = get_required_progs()?;
    
    let cli_args = Cli::parse();

    match cli_args.command {
        Commands::Build {
            fasta,
            gtf,
            rlen,
            output,
            spliced,
            unspliced,
            dedup,
            sparse,
            mut threads,
        } => {
            let r = run_fun!(mkdir -p $output)?;

            let info_file = output.join("run_info.json");
            let run_info = json!({
                "command" : "index",
                "version_info" : rp,
                "args" : {
                    "fasta" : fasta,
                    "gtf" : gtf,
                    "rlen" : rlen,
                    "output" : output,
                    "spliced" : spliced,
                    "unspliced" : unspliced,
                    "dedup" : dedup,
                    "sparse" : sparse,
                    "threads" : threads
                }
            });

            std::fs::write(&info_file, serde_json::to_string_pretty(&run_info).unwrap())
                .with_context(|| format!("could not write {}", info_file.display()))?;

            let outref = output.join("ref");
            let r = run_fun!(mkdir -p $outref)?;

            let mut cmd = std::process::Command::new(format!("{}", rp.pyroe.unwrap().display()));
            // we will run the make-splici command
            cmd.arg("make-splici");

            // if the user wants to dedup output sequences
            if dedup {
                cmd.arg(String::from("--dedup-seqs"));
            }
           
            // extra spliced sequence
            match spliced {
                Some(es) => {
                    cmd.arg(String::from("--extra-spliced"));
                    cmd.arg(format!("{}", es.display()));
                }
                None => {}
            }
           
            // extra unspliced sequence
            match unspliced {
                Some(eu) => {
                    cmd.arg(String::from("--extra-unspliced"));
                    cmd.arg(format!("{}", eu.display()));
                }
                None => {}
            }

            cmd.arg(fasta)
                .arg(gtf)
                .arg(format!("{}", rlen))
                .arg(&outref);
            let _cres = cmd.output()?;

            let mut salmon_index_cmd =
                std::process::Command::new(format!("{}", rp.salmon.unwrap().display()));
            let ref_prefix = format!("splici_fl{}.fa", rlen - 5);
            let ref_seq = outref.join(ref_prefix);

            let output_index_dir = output.join("index");
            salmon_index_cmd
                .arg("index")
                .arg("-i")
                .arg(output_index_dir)
                .arg("-t")
                .arg(ref_seq);

            // if the user requested a sparse index.
            if sparse {
                salmon_index_cmd.arg("--sparse");
            }

            // if the user requested more threads than can be used
            if let Ok(max_threads_usize) = std::thread::available_parallelism() {
                let max_threads = max_threads_usize.get() as u32;
                if threads > max_threads {
                    println!(
                        "The maximum available parallelism is {}, but {} threads were requested.",
                        max_threads, threads
                    );
                    println!("setting number of threads to {}", max_threads);
                    threads = max_threads;
                }
            }
           
            salmon_index_cmd
                .arg("--threads")
                .arg(format!("{}", threads));

            salmon_index_cmd.output().expect("failed to run salmon index");
        }
        Commands::Quant { index } => {
            println!("index is {}", index);
        }
    }
    Ok(())
}
