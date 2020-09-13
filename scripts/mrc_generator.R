
keepBSgenomeSequences <- function(genome, seqnames) {
  stopifnot(all(seqnames %in% seqnames(genome)))
  genome@user_seqnames <- setNames(seqnames, seqnames)
  genome@seqinfo <- genome@seqinfo[seqnames]
  genome
}

generate_mrcs <- function(df_sites_meta) {
  hg38 <- BSgenome.Hsapiens.UCSC.hg38
  sequences_to_keep <- paste0("chr", c(1:22, "X", "Y"))
  hg38 <- keepBSgenomeSequences(hg38, sequences_to_keep)  
  
  
  mrcs <- get_N_MRCs(df_sites_meta %>% dplyr::select(siteID, gender), hg38)
  mrcs <- merge(mrcs, df_sites_meta[,c("siteID", "sampleName", "refGenome")])
  mrcs$type <- "match"
  
  mrcs <- merge(mrcs,
                df_sites_meta[,c("siteID", "sampleName", "refGenome", "gender")]) %>%
    dplyr::select(-refGenome)
  
  
  sites <- df_sites_meta %>%
    dplyr::rename(chr = seqnames) %>% 
    dplyr::select(siteID, sampleName, gender, chr, strand, position, type)
  
  sites_mrcs <- rbind(sites, mrcs) %>%
    mutate(position = as.numeric(position))
  
  sites_mrcs <- makeGRanges(sites_mrcs,
                            soloStart = TRUE,
                            chromCol = 'chr',
                            strandCol = 'strand',
                            startCol = 'position')
  
  newSeqInfo <- seqinfo(hg38)
  seqInfo.new2old <- match(seqnames(newSeqInfo), seqnames(seqinfo(sites_mrcs)))
  seqinfo(sites_mrcs, new2old=seqInfo.new2old) <- newSeqInfo
  
  return(sites_mrcs)
}






