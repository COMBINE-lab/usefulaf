# Useful utilities for single-cell processing with alevin-fry

[Alevin-fry](https://github.com/COMBINE-lab/alevin-fry) is a fast, accurate and memory-frugal tool for preprocessing single-cell and single-nucleus RNA-seq data.  You can read more about alevin-fry in [its pre-print](https://www.biorxiv.org/content/10.1101/2021.06.29.450377v1).

This respoistory contains scripts, functions and utilities that are useful for preparing data for processing with alevin-fry, as well as for reading alevin-fry data into other packages for downstream analysis.

The different utilities are broken down in this repository by the language in which they are written (right now, Python, R and bash).  A brief listing of 
the available utilities currently in the repository is:

### R language 

* `make_splici_txome()` — A function to build a spliced + intron (_splici_) reference for indexing and quantification with `alevin-fry`. This function is available in the [`roe`](https://github.com/COMBINE-lab/roe) package, which can be installed by following [this instruction](https://github.com/COMBINE-lab/roe#installlation). 
* `emptyDropsCellRanger()` — An implementation of the hybrid UMI count filtering and [`emptyDrops`](https://github.com/MarioniLab/DropletUtils) used by CellRanger (and subsequently by [STARsolo](https://github.com/alexdobin/STAR)). This function is available in the [DropletUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html) BioConductor package.
* `loadFry()` — Contains a function to load `alevin-fry` output (including from USA mode quantification) into a [`SingleCellExperiment`](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) object. This function is available in the [fishpond](https://bioconductor.org/packages/release/bioc/html/fishpond.html) BioConductor package.

### Python language

* `load_fry.py` — Contains a Python function `load_fry` which is intended to load `alevin-fry` output (including from USA mode quantification) into a [`Scanpy`](https://github.com/theislab/scanpy) object.

### Bash

* `get_10x_permit_lists.sh` — Provides a script to download the 10x chromium v2 or v3 permit lists.
* `simpleaf` — Provides a script to run the entire `salmon -> alevin-fry (generate-permit-list > collate > quant)` pipeline, though providing only a simplified set of options.

-----------------

## Using simpleaf

  The `simpleaf` script that resides in the `bash` subdirectory is intended to simply the running of `alevin-fry` in common usage scenarios.  By limiting some of the different options that can be set, it provides a streamlined way to build the splici reference and index in a single command, as well as to process an experiment from raw FASTQ files to a count matrix in a single command.
  
  To work properly, `simpleaf` has a few requirements.  First, it should be run from *within* the `bash` subdirectory of this repository.  This is because it currently makes assumptions about the relative paths of the scripts `get_10x_permit_lists.sh` and `build_splici_ref.R`.  Additionally, the following environment variables are used within `simpleaf`:
  
   * `ALEVIN_FRY_HOME` **REQUIRED** — This directory will be used for persistent configuration and small file (<1G) storage between runs.  If you provide a directory and it doesn't exist, it will be created.  It is easiest to just set this in your enviornment globally so that the same home can be used over many runs without you having to provide the variable explicitly each time.  A good choice for this variable might be something like `~/.alevin_fry_home`.
   
   * `SALMON_BIN` **OPTIONAL** — This should provide the path to a `salmon` executable of version >= 1.5.1.  If not provided, the script will assume it can simply invoke `salmon` in the current enviornment.
   
   * `FRY_BIN` **OPTIONAL** — This should provide the path to a `alevin-fry` executable of version >= 0.4.0.  If not provided, the script will assume it can simply invoke `alevin-fry` in the current enviornment.
   
   * `TIME_BIN` **OPTIONAL** — This should provide the path to a [GNU time](https://www.gnu.org/software/time/) executable; this is different from the shell `time` command, and on most linux systems exists at `/usr/bin/time`.  If this variable is not provided, the script will assume it can use `/usr/bin/time`.  On OSX systems, you should install GNU time explicitly.  This can be done with [conda](https://anaconda.org/conda-forge/time) or homebrew.
  
  The `simpleaf` script has two sub-commands:
  
  * `index` — The `index` command will take a reference genome FASTA and GTF as input, build a splici reference using the `build_splici_ref.R` script, and then build a sparse `salmon` index on the resulting reference. **Note**: The `index` command requires the `Rscript` executable to be in the path, as well as all of theR packages that are required by `build_splici_ref.R`. The relevant options (which you can obtain by running `./simpleaf index -h`) are:
  
  ```{bash}
  Usage: ./simpleaf index [options]
        options:
         -f, --fasta REQUIRED genome reference FASTA file
         -g, --gtf REQUIRED GTF file with gene annotations
         -l, --rlen REQUIRED the target read length the index will be built for
         -o, --output REQUIRED path to output directory (will be created if it doesn't exist)
         -s, --spliced OPTIONAL path to FASTA file with extra spliced sequence to add to the index
         -u, --unspliced OPTIONAL path to FASTA file with extra unspliced sequence to add to the index
         -d, --dedup FLAG OPTIONAL deduplicate identical sequences inside the R script when building the splici reference
         -e, --dense FLAG OPTIONAL if this flag is passed, build the dense rather than sparse index for mapping
         -t, --threads OPTIONAL number of threads to use when running [default: min(16, num cores)]
         -h, --help display this help message
  ```
  
   * `quant` — The `quant` command takes as input the index, reads, and relevant information about the experiment (e.g. chemistry), and runs all of the steps of the `alevin-fry` pipeline, from mapping with `salmon` through `quant` with `alevin-fry`. The relevant options (which you can obtain by running `./simpleaf quant -h`) are:
  
  ```{bash}
  Usage: ./simpleaf quant [options]
        options:
         -1, --r1 REQUIRED comma separated list of left reads
         -2, --r2 REQUIRED comma separated list of right reads
         -i, --index REQUIRED path to a (sparse or dense) salmon splici index
         -o, --output REQUIRED path to output directory (will be created if it doesn't exist)
         -f, --fmode REQUIRED permit list filter mode, one of {knee, k, unfilt, u}
         -c, --chem REQUIRED chemistry of experiment, one of {v2, v3}
         -r, --res REQUIRED resolution strategy for alevin-fry, one of {cr-like, cr-like-em}
         -m, --t2g REQUIRED three-column txp-to-gene file to pass to alevin-fry quant command
         -t, --threads OPTIONAL number of threads to use when running [default: min(16, num cores)]
         -h, --help display this help message
  ```
