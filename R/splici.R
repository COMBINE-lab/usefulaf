make_splici_txome <- function(gtf_path, genome_path, read_length, flank_trim_length = 5, output_dir, dedup_seqs=FALSE) {
  # if you get some error from .get_cds_IDX, please try to rerun the code again
  # read length is the scRNA-seq read length
  # flank trim length is used to avoid marginal case when dealing with junctional reads 
  # assumes the following packages have been imported
  # eisaR, Biostrings, BSgenome, dplyr, stringr
  
  ########################################################################################################
  # Preprocessing
  ########################################################################################################
  
  suppressPackageStartupMessages({
    library(eisaR)
    library(Biostrings)
    library(BSgenome)
    library(stringr)
    library(GenomicFeatures)
  })
  
  if (!dir.exists(output_dir)) {
    dir.create(file.path(output_dir), showWarnings = FALSE)
  }  
  # make sure flank_length makes sense
  flank_length = read_length - flank_trim_length
  if (flank_length < 0 ){
    stop("flank trim length is larger than read length!")
  }
  # make sure gtf file exists
  if (!file.exists(gtf_path)) {
    stop("The following file does not exist: \n", gtf_path)
  }
  
  # make sure fasta file exist
  if (!file.exists(genome_path)) {
    stop("The following file does not exist: \n", genome_path)
  }
  file_name_prefix = paste0("transcriptome_splici_fl", flank_length)
  
  
  #########################################################################################################
  # Process gtf to get spliced and introns
  #########################################################################################################
  message("============processing gtf to get spliced and introns============")  
  # fl is the flank length, here we set it to 
  # the read length - 5 
  grl <- suppressWarnings(getFeatureRanges(
    gtf = file.path(gtf_path),
    featureType = c("spliced", "intron"), 
    intronType = "separate", 
    flankLength = flank_length, 
    joinOverlappingIntrons = TRUE, 
    verbose = TRUE
  ))
  
  #########################################################################################################
  # Get spliced related stuffs
  #########################################################################################################
  
  # spliced ranges has no dash in it
  spliced_grl = grl[sapply(strsplit(names(grl), "-"), length) == 1]
  
  #########################################################################################################
  # Get reduced introns
  #########################################################################################################
  
  # identify all introns and convert to GRanges
  intron_gr = unlist(grl[sapply(strsplit(names(grl), "-"), length) == 2])
  # group introns by gene, then collapse ovelaping ranges!
  intron_grl = reduce(split(intron_gr, intron_gr$gene_id))
  
  # clean txp names and gene names
  intron_gr <- BiocGenerics::unlist(intron_grl)
  intron_gr$exon_rank <- 1L
  intron_gr$transcript_id <- sapply(strsplit(names(intron_gr), "-"), function(x) x[1])
  intron_gr$gene_id <- intron_gr$transcript_id
  intron_gr$type <- "exon"
  intron_gr$transcript_id <- gsub(
    paste0("-I", "."), "-I", 
    make.unique(paste0(intron_gr$transcript_id, "-I")),
    fixed = TRUE
  )
  intron_gr$gene_id <- paste0(intron_gr$gene_id, "-I")
  intron_gr$exon_id <- intron_gr$transcript_id
  names(intron_gr) <- NULL
  mcols(intron_gr) <- 
    S4Vectors::mcols(intron_gr)[, c("exon_id", "exon_rank", 
                                    "transcript_id", "gene_id", "type")]
  # remake intron GRangesList
  intron_grl <- BiocGenerics::relist(intron_gr, lapply(
    structure(seq_along(intron_gr), 
              names = intron_gr$transcript_id), function(i) i))
  
  
  #########################################################################################################
  # extract sequences from genome
  #########################################################################################################
  
  message("============extracting spliced and intron sequences from genome============")  
  
  # load the genome sequence
  x <- Biostrings::readDNAStringSet(file.path(genome_path))
  # fix the names
  names(x) <- sapply(strsplit(names(x), " "), .subset, 1)
  
  
  grl = c(spliced_grl, intron_grl)
  
  # make sure introns don't out of boundary
  seqlevels(grl) <- seqlevels(x)
  seqlengths(grl) <- suppressWarnings(seqlengths(x)) 
  grl <- trim(grl)
  
  seqs <- GenomicFeatures::extractTranscriptSeqs(
    x = x,
    transcripts = grl
  )
  
  # If having duplicated sequences, only keep one
  if(dedup_seqs) {
    seqs = unique(seqs)
    grl = grl[names(seqs)]
  }
  
  # save some space
  rm(x)
  #########################################################################################################
  # process final outputs
  #########################################################################################################
  message("Writing outputs...")  
  
  df <- getTx2Gene(grl)
  write.table(df, file.path(output_dir, paste0(file_name_prefix, "_t2g.tsv")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
  df <- df %>%
    dplyr::mutate(gene_id = stringr::word(gene_id, 1, sep = '-'),
                  status = ifelse(stringr::str_detect(transcript_id, '-'), 'U', 'S'))
  
  writeXStringSet(seqs, file.path(output_dir, paste0(file_name_prefix, ".fa")), format = "fasta")
  write.table(df, file.path(output_dir, paste0(file_name_prefix, "_t2g_3col.tsv")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
  
  message("Done.")  
}
