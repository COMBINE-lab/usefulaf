# Usefulaf: An all-in-one Docker/Singularity image for single-cell processing with alevin-fry

[`Usefulaf`](https://hub.docker.com/r/combinelab/usefulaf/tags) is an all-in-one Docker/Singularity image for single-cell processing with [Alevin-fry](https://github.com/COMBINE-lab/alevin-fry)([paper](https://www.nature.com/articles/s41592-022-01408-3)). It includes the all tools you need to turn your FASTQ files into a count matrix and then load it into your favorite analysis environment. Specifically, this image includes:

- [`simpleaf`](https://github.com/COMBINE-lab/simpleaf): A simplified interface to indexing and quantifying with `alevin-fry`.
- [`pyroe`](https://github.com/COMBINE-lab/pyroe): An alevin-fry utility python package for building splici references, converting alevin-fry output formats, loading count matrix in Python, adding gene names (instead of just gene IDs) to output matrices, etc.
- [`fishpond::loadFry()`](https://rdrr.io/github/mikelove/fishpond/man/loadFry.html): A R function for loading count matrix as [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) object.

For processing data simply using the `usefulaf` image, check our latest tutorial [here](https://combine-lab.github.io/alevin-fry-tutorials/2021/quickstart-usefulaf-singularity/).

For pulling the Singularity image, please run the following code in bash. Note that the image is $\sim 1.65$ GB.

```bash
# if you use Docker
$ docker pull combinelab/usefulaf:latest

# if you use Singularity
$ singularity pull docker://combinelab/usefulaf:latest

```

## Usefulaf history

[Alevin-fry](https://github.com/COMBINE-lab/alevin-fry) is a fast, accurate, and memory-frugal tool for preprocessing single-cell and single-nucleus RNA-seq data. You can read more about alevin-fry in alevin-fry [pre-print](https://www.biorxiv.org/content/10.1101/2021.06.29.450377v2), and [paper](https://www.nature.com/articles/s41592-022-01408-3).

This repository was created initially with scripts, functions, and utilities that are useful for preparing data for processing with alevin-fry, as well as for reading alevin-fry data into other packages for downstream analysis. It also accompanies a Docker/Singularity container containing all of this relevant software in one place. However, as `alevin-fry` has continued to grow, all of that relevant functionality found its way into other, more stable and permanent homes (e.g. [`pyroe`](https://github.com/COMBINE-lab/pyroe) for splici reference construction and loading data in Python, [`roe`](https://github.com/COMBINE-lab/roe) for splici reference construction in R and [`fishpond`](https://bioconductor.org/packages/release/bioc/html/fishpond.html) for loading data in `R`). Finally, this repository also contained a bash script called `simpleaf` to simplify common workflows with `alevin-fry`. That, too, has evolved into its own (much more feature-rich and comprehensive) tool, living in its own repository ([`simpleaf`](https://github.com/COMBINE-lab/simpleaf)).

As such, all the scripts and functions in this repository have been retired. However, as usefulaf is still the only place that provides all these functionalities, we decided to turn [`usefulaf`] as an all-in-one [Docker/Singularity image](https://hub.docker.com/r/combinelab/usefulaf/tags) that makes use of all those new tools listed above. That has replaced the older `usefulaf` image that made use of the varied assortment of scripts and tools hosted in this repository. 
