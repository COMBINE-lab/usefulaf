

# assumes that SingleCellExperiment has already been loaded

load_fry <- function(frydir, which_counts = c('S','A'), verbose = FALSE) {
  library(rjson)
  library(Matrix)
  # read in metadata
  meta_info = fromJSON(file = file.path(frydir, "meta_info.json"))
  ng = meta_info$num_genes
  usa_mode = meta_info$usa_mode
  
  if(usa_mode) {
    if (length(which_counts) == 0){
      stop("Please at least provide one status in 'U' 'S' 'A' ")
    }
    if (verbose){
      message("processing input in USA mode, will return ", paste(which_counts, collapse = '+'))
    }
  } else if(verbose) {
    message("processing input in standard mode, will return spliced count")
  }

  # read in count matrix
  af_raw = readMM(file = file.path(frydir, "alevin", "quants_mat.mtx"))
  # if usa mode, each gene gets 3 rows, so ng/3
  if(usa_mode) {
    ng = as.integer(ng/3)
  }
  # read in gene name file and cell barcode file
  afg = read.csv(file.path(frydir, "alevin", "quants_mat_cols.txt"), strip.white = TRUE, header = FALSE, nrows = ng, col.names = c("gene_ids"))
  afc = read.csv(file.path(frydir, "alevin", "quants_mat_rows.txt"), strip.white = TRUE,header = FALSE,col.names = c("barcodes"))

  # if in usa_mode, sum up counts in different status according to which_counts
  if (usa_mode) {
    rd = list("S" = seq(1, ng), "U" =  seq(ng+1, 2*ng), "A" =  seq(2*ng+1, 3*ng))
    o = af_raw[, rd[[which_counts[1]]]]
    for (wc in which_counts[-1]) {
      o = o + af_raw[, rd[[wc]]]
    }
  } else {
    o = af_raw
  }
  
  # create SingleCellExperiment object
  sce <- SingleCellExperiment(list(counts = t(o)),
                              colData = afc,
                              rowData = afg,
  )
  sce
}