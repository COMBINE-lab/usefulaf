#!/usr/bin/env Rscript

# usage :
# $ ./build_splici_ref.R <path_to_genome_fasta> <path_to_gtf> <target_read_length> <output_dir>

# make sure the dependencies are installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install eisaR
BiocManager::install(c("eisaR","BSgenome"))

# install argparser
if (!requireNamespace("argparser", quietly = TRUE))
    install.packages("argparser")

# install roe from github
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

if (!requireNamespace("roe", quietly = TRUE))
    devtools::install_github("COMBINE-lab/roe")

# load packages
suppressPackageStartupMessages({
library(argparser)
library(roe)
})

# Create a parser
p <- arg_parser("Build a splici reference from a genome and GTF file.")

# Add command line arguments
# required arguments
p <- add_argument(p, "genome", help="path to genome.")
p <- add_argument(p, "gtf", help="GTF file with transcript annotations.")
p <- add_argument(p, "read-length", help="read length (determines flank size).",
	type="numeric")
p <- add_argument(p, "splici-dir", 
	help="output directory where splici reference information will be written.")

# optional arguments
p <- add_argument(p, "--flank-trim-length", 
	help="determines the amount subtracted from the read length to get the flank length.",
	type="numeric",
	default=5)
p <- add_argument(p, "--extra-spliced",
	help="path to FASTA file with any extra spliced sequences to include.")
p <- add_argument(p, "--extra-unspliced",
	help="path to FASTA file with any extra unspliced sequences to include.")
p <- add_argument(p, "--dedup-seqs", 
	help="deduplicate identical sequences before writing to output FASTA.", 
	flag=TRUE)

# Parse the command line arguments
argv <- parse_args(p)

# Set NAs to NULLs and call the function
if (is.na(argv$extra_spliced)) {
	argv$extra_spliced <- NULL
}
if (is.na(argv$extra_unspliced)) {
	argv$extra_unspliced <- NULL
}

make_splici_txome(gtf_path=argv$gtf,
		  genome_path=argv$genome,
		  read_length=argv$read_length,
		  flank_trim_length=argv$flank_trim_length,
		  output_dir=argv$splici_dir,
		  extra_spliced = argv$extra_spliced,
		  extra_unspliced = argv$extra_unspliced,
		  dedup_seqs = argv$dedup_seqs)
