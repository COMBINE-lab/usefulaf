# Useful utilities for single-cell processing with alevin-fry

[Alevin-fry](https://github.com/COMBINE-lab/alevin-fry) is a fast, accurate and memory-frugal tool for preprocessing single-cell and single-nucleus RNA-seq data.  You can read more about alevin-fry in [its pre-print](https://www.biorxiv.org/content/10.1101/2021.06.29.450377v1).

This respoistory contains scripts, functions and utilities that are useful for preparing data for processing with alevin-fry, as well as for reading alevin-fry data into other packages for downstream analysis.

The different utilities are broken down in this repository by the language in which they are written (right now, Python, R and bash).  A brief listing of 
the available utilities currently in the repository is:

### R language 

* `build_splici_ref.R` — A script to build a spliced + intron (splici) ref for indexing and quantification with `alevin-fry`.
* `splici.R` — Contains the `make_splici_txome` function, which is the function called by the `build_splici_ref.R` wrapper script.  If you want to build a splici reference programatically in R code, you can use this function.
* `cellRangerLikeEmptyDrops.R` — An implementation of the hybrid UMI count filtering and [`emptyDrops`](https://github.com/MarioniLab/DropletUtils) used by CellRanger (and subsequently by [STARsolo](https://github.com/alexdobin/STAR)). This R implementation is a translation of the implemntation in STARsolo, which itself was reverse-engineered from CellRanger. 
* `load_fry.R` — Contains a function to load `alevin-fry` output (including from USA mode quantification) into a [`SingleCellExperiment`](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) object.

### Python language

* `load_fry.py` — Contains a Python function `load_fry` which is intended to load `alevin-fry` output (including from USA mode quantification) into a [`Scanpy`](https://github.com/theislab/scanpy) object.

### Bash

* `get_10x_permit_lists.sh` — Provides a script to download the 10x chromium v2 or v3 permit lists.
* `simpleaf` — Provides a script to run the entire `salmon -> alevin-fry (generate-permit-list > collate > quant)` pipeline, though providing only a simplified set of options.
  * Notes about `simpleaf`. `simpleaf` requires you to set a `ALEVIN_FRY_HOME` variable in your path.  This directory will be used for persistent configuration and small file (<1G) storage between runs.  If you provide a directory and it doesn't exist, it will be created.  There are 2 sub-commands, `index` and `quant`.  The `index` command will take a reference genome FASTA and GTF as input, build a splici reference using the `build_splici_ref.R` script, and then build a sparse `salmon` index on the resulting reference.  The `quant` command takes as input the index, reads, and relevant information about the experiment (e.g. chemistry), and runs all of the steps of the `alevin-fry` pipeline, from mapping with `salmon` through `quant` with `alevin-fry`.  The `index` command requires the `Rscript` executable to be in the path, as well as all of theR packages that are required by `build_splici_ref.R`.  The `quant` command requires that the `salmon` and `alevin-fry` execuitables be in the path OR that the variables `SALMON_BIN` and `FRY_BIN` exist in the environment.  The `simpleaf` script will check the versions of `salmon` and `alevin-fry`, but recent versions of both are required (`salmon` v >= 1.5.1, `alevin-fry` v >= 0.4.0). 
