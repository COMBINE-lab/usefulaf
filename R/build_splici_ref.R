#!/usr/bin/env Rscript

# usage :
# $ ./build_splici_ref.R <path_to_genome_fasta> <path_to_gtf> <target_read_length> <output_dir>

# install BioC depedencies if necessary
if ( (!requireNamespace("eisaR", quietly = TRUE)) || 
     (!requireNamespace("BSgenome", quietly = TRUE)) || 
     (!requireNamespace("fishpond", quietly = TRUE)) ) {

  # install BioC itself, if we don't have it
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  BiocManager::install(c("eisaR","BSgenome","fishpond"))
}

# install argparser
if (!requireNamespace("argparser", quietly = TRUE))
    install.packages("argparser")

# install devtools 
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# install roe from github
if (!requireNamespace("roe", quietly = TRUE))
    devtools::install_github("COMBINE-lab/roe",dep = FALSE)

# load packages
suppressPackageStartupMessages({
library(argparser)
library(roe)
})

# Create a parser
p <- arg_parser("Build a splici reference from a genome and GTF file.")

# Add command line arguments
# required arguments
p <- add_argument(p, "genome", help="The path to a genome FASTA file.")
p <- add_argument(p, "gtf", help="The path to a gtf file.")
p <- add_argument(p, "read-length", help="The read length of the single-cell experiment being processed (determines flank size).",
	type="numeric")
p <- add_argument(p, "output-dir", 
	help="The output directory where splici reference files will be written.")

# optional arguments
p <- add_argument(p, "--flank-trim-length", 
	help="Determines the amount subtracted from the read length to get the flank length.",
	type="numeric",
	default=5)
p <- add_argument(p, "--filename-prefix", 
	help="The file name prefix of the generated output files.",
	default="splici")
p <- add_argument(p, "--extra-spliced",
	help="The path to an extra spliced sequence fasta file.")
p <- add_argument(p, "--extra-unspliced",
	help="The path to an extra unspliced sequence fasta file.")
p <- add_argument(p, "--dedup-seqs",
	help="A flag indicates whether identical sequences will be deduplicated.", 
	flag=TRUE)
p <- add_argument(p, "--no-flanking-merge", 
	help="A flag indicates whether flank lengths will be considered when merging introns.", 
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

make_splici_txome(gtf_path = argv$gtf,
				genome_path = argv$genome,
				read_length = argv$read_length,
				output_dir = argv$output_dir,
				flank_trim_length = argv$flank_trim_length,
				filename_prefix = argv$filename_prefix,
				extra_spliced = argv$extra_spliced,
				extra_unspliced = argv$extra_unspliced,
				dedup_seqs = argv$dedup_seqs,
				no_flanking_merge = argv$no_flanking_merge)
