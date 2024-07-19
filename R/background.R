# Background functions for QTLretrievR that will be called by the primary functions

## usethis namespace: start
#' @importFrom magrittr %>%
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom DT datatable
#' @import qtl2
#' @importFrom intermediate mediation.scan
## usethis namespace: end

split_map <- function(map, chr_names = NULL){
  map <- reorder_map_table(map, chr_names = chr_names)
  pos <- as.numeric(map[, 2])
  chr <- map[, 1]
  uchr <- unique(chr)
  names(pos) <- rownames(map)
  lapply(split(pos, factor(chr, uchr)), sort)
}

reorder_map_table <- function(map_tab, chr_col = 1, pos_col = 2, chr_names = NULL){
  chr <- map_tab[, chr_col]
  if (is.null(chr_names))
    chr_names <- unique(chr)
  chr <- factor(chr, levels = chr_names)
  pos <- suppressWarnings(as.numeric(map_tab[, pos_col]))
  map_tab[order(chr, pos, seq_len(nrow(map_tab))), , drop = FALSE]
}

rankZ <- function(x){
  x <- rank(x, na.last = "keep", ties.method = "average")/(sum(!is.na(x)) + 1)
  qnorm(x)
}

probs_3d_to_qtl2 <- function(probs){
  uchroms <- c(as.character(1:19), "X")
  # Convert to qtl2 genoprobs format
  # Similar to qtl2convert::probs_doqtl_to_qtl2()
  markers <- dimnames(probs)[[3]]
  chroms <- sapply(strsplit(markers, "_"), "[[", 1)
  newprobs <- vector("list", length(uchroms))
  names(newprobs) <- uchroms
  for (chrom in uchroms) newprobs[[chrom]] <- probs[, , chroms == chrom]
  attr(newprobs, "is_x_chr") <- c(rep(FALSE, length(uchroms)-1),TRUE)
  attr(newprobs, "crosstype") <- "DO"
  attr(newprobs, "alleles") <- c("A","B","C","D","E","F","G","H")
  attr(newprobs, "alleleprobs") <- TRUE
  class(newprobs) <- c("calc_genoprob", "list")
  newprobs
}

batchmap <- function(nbatch, exprZ, kinship_loco, genoprobs, covar, tissue, gmap, thrA=5, thrX=5, n.cores = 4, ...){
  # Do mapping in batches. I will only save significant peaks, for efficiency.
  nn <- ncol(exprZ)
  ss <- round(seq(0, nn, length.out=nbatch + 1))
  peaks <- list()
  xx <- parallel::detectCores()
  cl <- parallel::makeCluster(floor(xx/n.cores))
  doParallel::registerDoParallel(cl)
  foreach::foreach(i = 1:nbatch) %dopar% {
    start <- ss[i] + 1
    end <- ss[i + 1]
    cat(sprintf("batch %d: %d-%d\n", i, start, end))
    out <- qtl2::scan1(genoprobs, exprZ[, start:end, drop=FALSE],
                       kinship_loco, addcovar=covar[,,drop=FALSE], cores=n.cores, ...)
    peaks[[i]] <- qtl2::find_peaks(out, gmap, drop=1.5,
                                   threshold=thrA, thresholdX=thrX)   # returns a long & tidy dataset
  }
  parallel::stopCluster(cl)
  do.call('rbind', peaks) %>% dplyr::select(-lodindex) %>%
    dplyr::rename(phenotype=lodcolumn, peak_chr=chr, peak_cM=pos)
}

interp_bp <- function(df) {
  chroms <- c(as.character(1:19), "X")
  df <- arrange(df, peak_chr, peak_cM)
  peak_gpos <- select(df, peak_chr, peak_cM)
  chr <- peak_gpos$peak_chr
  f <- factor(chr, chroms)
  peak_gcoord_list <- split(peak_gpos$peak_cM, f)
  peak_pcoord_list <- qtl2::interp_map(peak_gcoord_list, gmap, pmap)
  df$interp_bp_peak <- unsplit(peak_pcoord_list, f)
  df
}

create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                rownames = FALSE,
                filter="top",
                options = list(dom = 'Blfrtip',
                               buttons = c('copy', 'csv', 'excel'),
                               pageLength = 5,
                               scrollX= TRUE
                ))

}

standardize <- function(expr, exprZ, details, tissues = c()){
  tissue_samp <- list()
  for(tissue in tissues){
    rownames(expr[[tissue]]) <- gsub("_.*","", rownames(expr[[tissue]]))
    colnames(exprZ[[tissue]]) <- gsub("_.*","", colnames(exprZ[[tissue]]))
    tissue_samp[[tissue]] <- sample_details[which(sample_details$ID %in% colnames(expr[[tissue]])),]
  }
  return(list(expr, exprZ, tissue_samp))
}

batchmediate <- function( n, z_thres = -2,  pos_thres = 10, QTL.peaks, med_annot, QTL.mediator, targ_covar, QTL.target, probs, ...){

  med.scan <- list()

  #for( l in 1:nbatch){

  start <- ss[n]+1
  end   <- ss[n+1]
  lod.peaks <- QTL.peaks[start:end,]
  cat(sprintf("batch %d: %d-%d\n", n, start, end))

  for(i in 1:nrow(lod.peaks) ){

    marker    <- map_dat2 %>%
      dplyr::mutate(pos=as.numeric(pos_bp)) %>%
      dplyr::filter( abs(pos - lod.peaks$peak_cM[i]) == min(abs(pos - lod.peaks$peak_cM[i])))
    qtl.chr   <- marker$chr
    qtl.pos   <- marker$pos_bp/1e06
    annot     <- med_annot %>% mutate( middle_point = pos)
    geno      <- qtl2::pull_genoprobpos(probs,marker$marker)
    geno      <- geno[rownames(geno) %in% rownames(QTL.target),]
    target    <- lod.peaks$phenotype[i]

    med <- intermediate::mediation.scan(target     = QTL.target[, target, drop = FALSE],
                                        mediator   = QTL.mediator,
                                        annotation = annot,
                                        qtl.geno   = geno,
                                        covar      = targ_covar,
                                        verbose    = FALSE,
                                        method     = "double-lod-diff")

    med <- med %>%
      mutate(
        target_id = lod.peaks$phenotype[i],
        qtl_lod = lod.peaks$lod[i],
        qtl_chr = lod.peaks$peak_chr[i]
      ) %>%
      select(
        target_id,
        qtl_lod,
        qtl_chr,
        mediator = symbol,
        mediator_id ,
        mediator_chr = chr,
        mediator_midpoint = middle_point,
        LOD)



    med.scan[[i]] <- med

    print(i)
  }

  return(med.scan)
  #}
}

subset_probs <- function(this_probs, this_chrom, this_markers) {
  att <- attributes(this_probs)
  att$names <- this_chrom
  att$is_x_chr <- setNames(FALSE, this_chrom)
  assert_that(all(this_markers %in% dimnames(this_probs[[this_chrom]])[[3]]))
  newprobs <- list(this_probs[[this_chrom]][, , this_markers, drop=FALSE])
  names(newprobs) <- this_chrom
  attributes(newprobs) <- att
  newprobs
}

`%notin%` <- Negate(`%in%`)
