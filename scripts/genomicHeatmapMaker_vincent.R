packs <- c("colorspace", "hiAnnotator", "GCcontent", "BSgenome")
null <- suppressMessages(sapply(packs, library, character.only = TRUE))

getRefSeq_genes <- function(reference_genome) {
  # Identify table and track names
  
  
  bsession <- makeUCSCsession(reference_genome)
  trk_names <- names(trackNames(ucscTableQuery(bsession)))
  tbl_names <- tableNames(ucscTableQuery(bsession))
  
  trk <- grep("RefSeq", trk_names, value = TRUE)
  trk <- trk[trk %in% c("NCBI RefSeq", "RefSeq Genes")]
  stopifnot(length(trk) == 1)
  
  tbl <- grep("refGene", tbl_names, value = TRUE)
  stopifnot(length(tbl) == 1)
  
  refSeq <- makeGRanges(
    getUCSCtable(tbl, trk, freeze=reference_genome)
    # freeze=reference_genome
  )
  
  return(refSeq)
}

get_cpg <- function(referenceGenome, reference_genome_sequence) {
  #CpG_islands <- getCpG_islands(referenceGenome)
  CpG_data <- cpg <- getUCSCtable(
    "cpgIslandExt", "CpG Islands", freeze = referenceGenome
  )
  
  CpG_data <- CpG_data[
    CpG_data$chrom %in% 
      seqnames(GenomeInfoDb::seqinfo(reference_genome_sequence)),
  ]
  
  CpG_islands <- GenomicRanges::GRanges(
    seqnames = CpG_data$chrom,
    ranges = IRanges::IRanges(
      start = CpG_data$chromStart, end = CpG_data$chromEnd
      ),
    strand = "*",
    seqinfo = GenomeInfoDb::seqinfo(reference_genome_sequence)
  )
  
  mcols(CpG_islands) <- CpG_data  
  
  return(CpG_islands)
}

get_oncogene_from_file <- function(filename) {
  onco <- read.csv(filename, header=FALSE, stringsAsFactors=FALSE)
  as.character(onco$V1)
}

getPositionalValuesOfFeature <- function(sites, genomicData) {
  #### Boundary Distances #### Nirav Malani code TODO: refactor into several functions
  ## (refSeq boundary.dist), Start (refSeq start.dist), non-width (), General (general.width)
  ## when inGene is FALSE then set following: ref.left.pos, ref.right.pos, ref.left.strand, ref.right.strand
  ## when inGene is TRUE then set following: ref.start.pos, ref.end.pos, ref.gene.strand
  
  ## prepare the new columns ##
  colnam <- paste("ref", c("left.pos", "right.pos", "left.strand", "right.strand", 
                           "start.pos", "end.pos", "gene.strand"), sep=".") 
  mcols(sites)[colnam] <- NA
  
  ## add the respective columns as needed ##
  ## beware: precede returns range which is following the query and
  ## follow returns the range which is preceding the query!
  ## so do a switcheroo in terms of extracting the start & stop ##
  left <- follow(sites, genomicData, ignore.strand=TRUE)
  left[is.na(left) | sites$within_refSeq_gene] <- NA
  rows <- na.omit(left)
  sites$ref.left.pos[!is.na(left)] <- end(genomicData[rows])
  sites$ref.left.strand[!is.na(left)] <- as.character(strand(genomicData[rows]))
  
  right <- precede(sites, genomicData, ignore.strand=TRUE)
  right[is.na(right) | sites$within_refSeq_gene] <- NA
  rows <- na.omit(right)
  sites$ref.right.pos[!is.na(right)] <- start(genomicData[rows])
  sites$ref.right.strand[!is.na(right)] <- as.character(strand(genomicData[rows]))
  
  inIt <- findOverlaps(sites, genomicData, ignore.strand=TRUE, select="arbitrary")
  inIt[is.na(inIt) | !sites$within_refSeq_gene] <- NA
  rows <- na.omit(inIt)
  sites$ref.start.pos[!is.na(inIt)] <- start(genomicData[rows])
  sites$ref.end.pos[!is.na(inIt)] <- end(genomicData[rows])
  sites$ref.gene.strand[!is.na(inIt)] <- as.character(strand(genomicData[rows]))
  
  sites$boundary.dist <-
    eval(expression(pmin((ref.end.pos-position)/(ref.end.pos-ref.start.pos),
                         (position-ref.start.pos)/(ref.end.pos-ref.start.pos),
                         (ref.right.pos-position)/(ref.right.pos-ref.left.pos),
                         (position-ref.left.pos)/(ref.right.pos-ref.left.pos),
                         na.rm=T)), mcols(sites))
  
  sites$start.dist <-
    eval(expression(pmin(ifelse(ref.gene.strand=="-",
                                (ref.end.pos-position)/(ref.end.pos-ref.start.pos),
                                (position-ref.start.pos)/(ref.end.pos-ref.start.pos)),
                         ifelse(ref.right.strand=="-",
                                (ref.right.pos-position)/(ref.right.pos-ref.left.pos),
                                NA),
                         ifelse(ref.left.strand=="+",
                                (position-ref.left.pos)/(ref.right.pos-ref.left.pos),
                                NA),na.rm=T)), mcols(sites))
  
  sites$general.width <- eval(expression(pmin(ref.end.pos-ref.start.pos, 
                                              ref.right.pos-ref.left.pos,na.rm=T)),
                              mcols(sites))
  sites$gene.width <- eval(expression(ref.end.pos-ref.start.pos ), mcols(sites))
  
  meta <- mcols(sites)
  meta <- meta[ , ! (names(meta) %in% colnam)]
  mcols(sites) <- meta
  
  sites 
}

from_counts_to_density <- function(sites, column_prefix, window_size) {
  metadata <- mcols(sites)
  sapply(seq(window_size), function(i) {
    val <- window_size[i]
    name <- names(window_size)[i]
    column_name <- paste0(column_prefix, ".", name)
    metadata[[column_name]] <<- metadata[[column_name]]/val
  })
  mcols(sites) <- metadata
  sites
}

get_annotation_columns <- function(sites) {
  granges_column_names <- c("seqnames", "start", "end", "width", "strand")
  int_site_column_names <- c("siteID", "sampleName", "chr", "strand", "position")
  required_columns <- unique(c(
    granges_column_names, int_site_column_names, "type"))
  stopifnot(all(required_columns %in% names(sites)))
  setdiff(names(sites), required_columns)
}

# modified as.numeric
na.median <- function(x) {
    if (!is.matrix(x)) x <- as.matrix(x)
    na.count <- colSums(is.na(x))
    if (any(na.count != 0 ))
    {
      for (i in 1:ncol(x))
        if (na.count[i]>0){
          med <- median(as.numeric(x[,i]), na.rm=TRUE)
          x[is.na(x[,i]),i] <- med
        }
    }
    x
}

#modify gtsp to sampleName
sites_to_ROC_ordinary <- function(sites_mrcs, sampleName_GTSP, output_dir) {
  sites_mrcs <- as.data.frame(sites_mrcs)
  
  #write.table(sites_mrcs, file='sites_mrcs.gen')
  
  annotation_columns <- get_annotation_columns(sites_mrcs)
  
  #i <- sites_mrcs[c("type", annotation_columns, "sampleName")] 
  #message('sites_mrcs data frame:')
  #write.table(i)
  
  roc.res <- ROC.ORC(
    sites_mrcs[,"type"],
    na.median(sites_mrcs[,annotation_columns]),
    sites_mrcs[,"sampleName"],
    origin.levels = unique(as.character(sites_mrcs[,"sampleName"])))
  
  ROCSVG(roc.res, output_dir)
  saveRDS(roc.res, file = file.path(output_dir, "roc.res.rds"))
}