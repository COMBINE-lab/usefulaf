use anyhow::{anyhow, bail, Result};
use clap::{Parser, Subcommand};
use cmd_lib::run_fun;
use semver::{Version, VersionReq};
use std::env;
use std::ffi::OsStr;
use std::fs;
use std::path::PathBuf;
use std::process::Command;
use which::which;

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

// Holds the paths to the
// programs we'll need to run
// the tool.
struct ReqProgs {
    salmon: Option<PathBuf>,
    alevin_fry: Option<PathBuf>,
    pyroe: Option<PathBuf>,
}

fn search_for_executable(env_key: &str, prog_name: &str) -> Result<PathBuf> {
    match env::var(env_key) {
        Ok(p) => {
            return Ok(PathBuf::from(p));
        }
        Err(e) => {
            eprintln!("${} is unset {}, trying default path.", env_key, e);
            eprintln!(
                "If a satisfactory version is not found, consider setting the ${} variable.",
                env_key
            );
            match which(prog_name) {
                Ok(p) => {
                    println!("found `{}` in the PATH at {}", prog_name, p.display());
                    return Ok(p);
                }
                Err(e) => {
                    return Err(anyhow!(
                        "could not find `{}` in your path: {}",
                        prog_name,
                        e
                    ));
                }
            }
        }
    }
}

fn check_version_constraints<S1: AsRef<str>>(
    req_string: S1,
    prog_output: std::result::Result<String, std::io::Error>,
) -> Result<Version> {
    match prog_output {
        Ok(vs) => {
            let x = vs.split_whitespace();
            if let Some(version) = x.last() {
                let parsed_version = Version::parse(version).unwrap();
                let req = VersionReq::parse(req_string.as_ref()).unwrap();
                if req.matches(&parsed_version) {
                    return Ok(parsed_version);
                } else {
                    return Err(anyhow!(
                        "parsed version {:?} does not satisfy constraints {:?}",
                        version,
                        req
                    ));
                }
            }
        }
        Err(e) => {
            eprintln!("Error running salmon {}", e);
            return Err(anyhow!("could not parse program output"));
        }
    }
    Err(anyhow!("invalid version string"))
}

fn get_required_progs() -> Result<ReqProgs> {
    let mut rp = ReqProgs {
        salmon: None,
        alevin_fry: None,
        pyroe: None,
    };

    // First look for any environment variables
    // then check the path.
    rp.salmon = Some(search_for_executable("SALMON", "salmon")?);
    rp.alevin_fry = Some(search_for_executable("ALEVIN_FRY", "alevin-fry")?);
    rp.pyroe = Some(search_for_executable("PYROE", "pyroe")?);

    if let Some(salmon) = rp.salmon.clone() {
        let st = salmon.display().to_string();
        println!("st = {}", st);
        let sr = run_fun!($st --version);
        let v = check_version_constraints(">=1.5.1, <2.0.0", sr)?;
        println!("{:?}, {}", v, v);
    }

    if let Some(af) = rp.alevin_fry.clone() {
        let st = af.display().to_string();
        println!("st = {}", st);
        let sr = run_fun!($st --version);
        let v = check_version_constraints(">=0.4.1, <1.0.0", sr)?;
        println!("** {:?}, {}", v, v);
    }

    if let Some(pr) = rp.pyroe.clone() {
        let st = pr.display().to_string();
        println!("st = {}", st);
        let sr = run_fun!($st --version);
        let v = check_version_constraints(">=0.6.2, <1.0.0", sr)?;
        println!("** {:?}, {}", v, v);
    }
    Ok(rp)
}

fn main() -> anyhow::Result<()> {
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
            threads,
        } => {
            let outref = output.join("ref");
            let r = run_fun!(mkdir -p $outref)?;

            let mut cmd_str = format!("{} make-splici", rp.pyroe.unwrap().display());

            //let fp = format!(" --filename-prefix splici");
            //cmd_str += &fp;
            let dedup_str = if dedup {
                Some(String::from("--dedup-seqs"))
            } else {
                None
            };

            let spliced_str = match spliced {
                Some(es) => Some(format!("--extra-spliced {}", es.display())),
                None => None,
            };

            let unspliced_str = match unspliced {
                Some(eu) => Some(format!("--extra-unspliced {}", eu.display())),
                None => None,
            };

            let build_res = match (dedup_str, spliced_str, unspliced_str) {
                (Some(d), Some(s), Some(u)) => {
                    run_fun!(pyroe make-splici $d $s $u $fasta $gtf $rlen $outref)?
                }
                (Some(d), Some(s), None) => {
                    run_fun!(pyroe make-splici $d $s $fasta $gtf $rlen $outref)?
                }
                (Some(d), None, None) => run_fun!(pyroe make-splici $d $fasta $gtf $rlen $outref)?,
                (Some(d), None, Some(u)) => {
                    run_fun!(pyroe make-splici $d $u $fasta $gtf $rlen $outref)?
                }
                (None, Some(s), Some(u)) => {
                    run_fun!(pyroe make-splici $s $u $fasta $gtf $rlen $outref)?
                }
                (None, Some(s), None) => run_fun!(pyroe make-splici $s $fasta $gtf $rlen $outref)?,
                (None, None, Some(u)) => run_fun!(pyroe make-splici $u $fasta $gtf $rlen $outref)?,
                (None, None, None) => run_fun!(pyroe make-splici $fasta $gtf $rlen $outref)?,
            };

            println!("{:?}", build_res);
        }
        Commands::Quant { index } => {
            println!("index is {}", index);
        }
    }
    Ok(())
}
