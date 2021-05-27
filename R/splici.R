make_splici_txome <- function(gtf_path, genome_path, flank_length, output_dir) {
    # assumes the following packages have been imported
    # eisaR, Biostrings, BSgenome, dplyr, stringr
    dir.create(output_dir)

    gtf <- file.path(gtf_path)
    # fl is the flank length, here we set it to 
    # the read length - 5 (151 - 5bp = 146bp) 
    grl <- getFeatureRanges(
      gtf = gtf,
      featureType = c("spliced", "intron"), 
      intronType = "separate", 
      flankLength = flank_length, 
      joinOverlappingIntrons = FALSE, 
      verbose = TRUE
    )

    # load the genome sequence
    x <- Biostrings::readDNAStringSet(file.path(genome_path))
    # fix the names
    names(x) <- stringr::word(names(x), 1)

    seqlevels(grl) <- seqlevels(x)
    seqlengths(grl) <- seqlengths(x)

    grl <- trim(grl)

    seqs <- GenomicFeatures::extractTranscriptSeqs(
      x = x,
      transcripts = grl
    )

    seqs <- unique(seqs)
    grl <- grl[names(seqs)]

    df <- getTx2Gene(grl)
    write.table(df, file.path(output_dir, "t2g.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
    df <- df %>%
      dplyr::mutate(gene_id = stringr::word(gene_id, 1, sep = '-'),
                    status = ifelse(stringr::str_detect(transcript_id, '-'), 'U', 'S'))

    writeXStringSet(seqs, file.path(output_dir, "transcriptome_splici.fa"), format = "fasta")
    write.table(df, file.path(output_dir, "t2g_3col.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
}
