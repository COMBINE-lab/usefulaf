# Useful utilities for single-cell processing with alevin-fry

[Alevin-fry](https://github.com/COMBINE-lab/alevin-fry) is a fast, accurate and memory-frugal tool for preprocessing single-cell and single-nucleus RNA-seq data.  You can read more about alevin-fry in alevin-fry [pre-print](https://www.biorxiv.org/content/10.1101/2021.06.29.450377v2), and [paper](https://www.nature.com/articles/s41592-022-01408-3).

This respoistory was created initially with scripts, functions and utilities that are useful for preparing data for processing with alevin-fry, as well as for reading alevin-fry data into other packages for downstream analysis. It also accompanies a Docker/Singularity container containing all of this relevant software in one place.  However, as `alevin-fry` has continued to grow, all of that relevant functionality found it's way into other, more stable and permanent homes (e.g. [`pyroe`](https://github.com/COMBINE-lab/pyroe) for splici reference construction and loading data in python, and [`fishpond`](https://bioconductor.org/packages/release/bioc/html/fishpond.html) for loading data in `R`).  Finally, this repository also contained a bash script, called `simpleaf` to simplify common workflows with `alevin-fry`.  That, too, has evolved into it's own (much more feature-rich and comprehensive) tool, living in it's own repository ([`simpleaf`](https://github.com/COMBINE-lab/simpleaf)).

**As such** this repository has been retired.  If you are interested in the functionality that it provides, please look at the following resources instead:

* [`simpleaf`](https://github.com/COMBINE-lab/simpleaf) — for a simplified interface to indexing and quantifying with `alevin-fry`.

* [`pyroe`](https://github.com/COMBINE-lab/pyroe) (and [on pip](https://pypi.org/project/pyroe/)) — for building splici references, converting `alevin-fry` output formats, loading data in `alevin-fry`, adding gene names (instead of just gene IDs) to output matrices, etc.

* [`fishpond`](https://bioconductor.org/packages/release/bioc/html/fishpond.html) — for ingesting `alevin-fry` data in R to interact with R's single-cell analysis ecosystem.

Finally, we provide an all-in-one Docker/Singularity image called [`usefulaf`](https://hub.docker.com/r/combinelab/usefulaf/tags) that makes use of the newer `simpleaf` tool. That has replaced the older `usefulaf` image that made use of the varied assortment of scripts and tools hosted in this repository.  The updated tutorial describing how to process data simply using the `usefulaf` image can be found [here](https://combine-lab.github.io/alevin-fry-tutorials/2021/quickstart-usefulaf-singularity/).
