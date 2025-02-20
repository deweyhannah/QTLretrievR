#' Background functions for QTLretrievR that will be called by the primary functions
#' @importFrom magrittr %>%
#' @importFrom foreach foreach %dopar%
#' @importFrom dplyr rename select arrange
#' @importFrom DT datatable
#' @importFrom doParallel registerDoParallel
#' @importFrom stats setNames
#' @import qtl2
#' @importFrom intermediate mediation.scan
#' @importFrom stats qnorm
#'

split_map <- function(map, chr_names = NULL) {
  map <- reorder_map_table(map, chr_names = chr_names)
  pos <- as.numeric(map[, 2])
  chr <- map[, 1]
  uchr <- unique(chr)
  names(pos) <- rownames(map)
  lapply(split(pos, factor(chr, uchr)), sort)
}

reorder_map_table <- function(map_tab, chr_col = 1, pos_col = 2, chr_names = NULL) {
  chr <- map_tab[, chr_col]
  if (is.null(chr_names)) {
    chr_names <- unique(chr)
  }
  chr <- factor(chr, levels = chr_names)
  pos <- suppressWarnings(as.numeric(map_tab[, pos_col]))
  map_tab[order(chr, pos, seq_len(nrow(map_tab))), , drop = FALSE]
}

rankZ <- function(x) {
  x <- rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  qnorm(x)
}

probs_3d_to_qtl2 <- function(probs) {
  uchroms <- c(as.character(1:19), "X")
  # Convert to qtl2 genoprobs format
  # Similar to qtl2convert::probs_doqtl_to_qtl2()
  markers <- dimnames(probs)[[3]]
  chroms <- sapply(strsplit(markers, "_"), "[[", 1)
  newprobs <- vector("list", length(uchroms))
  names(newprobs) <- uchroms
  for (chrom in uchroms) newprobs[[chrom]] <- probs[, , chroms == chrom]
  attr(newprobs, "is_x_chr") <- c(rep(FALSE, length(uchroms) - 1), TRUE)
  attr(newprobs, "crosstype") <- "DO"
  attr(newprobs, "alleles") <- c("A", "B", "C", "D", "E", "F", "G", "H")
  attr(newprobs, "alleleprobs") <- TRUE
  class(newprobs) <- c("calc_genoprob", "list")
  newprobs
}

interp_bp <- function(df, genmap, physmap) {
  chroms <- c(as.character(1:19), "X")
  df <- dplyr::arrange(df, peak_chr, peak_cM)
  peak_gpos <- select(df, peak_chr, peak_cM)
  chr <- peak_gpos$peak_chr
  f <- factor(chr, chroms)
  peak_gcoord_list <- split(peak_gpos$peak_cM, f)
  peak_pcoord_list <- qtl2::interp_map(peak_gcoord_list, genmap, physmap)
  df$interp_bp_peak <- unsplit(peak_pcoord_list, f)
  df
}

create_dt <- function(x) {
  DT::datatable(x,
    extensions = "Buttons",
    rownames = FALSE,
    filter = "top",
    options = list(
      dom = "Blfrtip",
      buttons = c("copy", "csv", "excel"),
      pageLength = 5,
      scrollX = TRUE
    )
  )
}

standardize <- function(expr, exprZ, details, tissues = c()) {
  tissue_samp <- list()
  for (tissue in tissues) {
    rownames(expr[[tissue]]) <- gsub("_.*", "", rownames(expr[[tissue]]))
    colnames(exprZ[[tissue]]) <- gsub("_.*", "", colnames(exprZ[[tissue]]))
    tissue_samp[[tissue]] <- details[which(details$ID %in% colnames(expr[[tissue]])) ,] # & details$tissue == tissue), ]
  }
  std <- tibble::lst(expr, exprZ, tissue_samp)
  # names(std) <- c("expr","exprZ","tissue_samp")
  return(std)
}

batchmediate <- function(n, z_thres = -2, pos_thres = 10, QTL.peaks, med_annot, QTL.mediator, targ_covar, QTL.target, probs, ...) {
  med.scan <- list()

  start <- ss[n] + 1
  end <- ss[n + 1]
  lod.peaks <- QTL.peaks[start:end, ]
  cat(sprintf("batch %d: %d-%d\n", n, start, end))

  for (i in 1:nrow(lod.peaks)) {
    marker <- map_dat2 %>%
      dplyr::mutate(pos = as.numeric(pos_bp)) %>%
      dplyr::filter(abs(pos - lod.peaks$peak_cM[i]) == min(abs(pos - lod.peaks$peak_cM[i])))
    qtl.chr <- marker$chr
    qtl.pos <- marker$pos_bp / 1e06
    annot <- med_annot %>% mutate(middle_point = pos)
    geno <- qtl2::pull_genoprobpos(probs, marker$marker)
    geno <- geno[rownames(geno) %in% rownames(QTL.target), ]
    target <- lod.peaks$phenotype[i]

    med <- intermediate::mediation.scan(
      target = QTL.target[, target, drop = FALSE],
      mediator = QTL.mediator,
      annotation = annot,
      qtl.geno = geno,
      covar = targ_covar,
      verbose = FALSE,
      method = "double-lod-diff"
    )

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
        mediator_id,
        mediator_chr = chr,
        mediator_midpoint = middle_point,
        LOD
      )

    med.scan[[i]] <- med

    print(i)
  }

  return(med.scan)
}

subset_probs <- function(this_probs, this_chrom, this_markers) {
  att <- attributes(this_probs)
  att$names <- this_chrom
  att$is_x_chr <- setNames(FALSE, this_chrom)
  # assert_that(all(this_markers %in% dimnames(this_probs[[this_chrom]])[[3]]))
  newprobs <- list(this_probs[[this_chrom]][, , this_markers, drop = FALSE])
  names(newprobs) <- this_chrom
  attributes(newprobs) <- att
  newprobs
}

`%notin%` <- Negate(`%in%`)

peak_fun <- function(i,ss, exprZ, kinship_loco, genoprobs, covar, tissue, gmap, thrA = 5, thrX = 5, n.cores = 4) {
  # message("Started mapping")
  # message(timestamp())
  start <- ss[i] + 1
  end <- ss[i + 1]
  out <- qtl2::scan1(
    genoprobs,
    exprZ[,start:end,drop=F],
    kinship_loco,
    addcovar = covar[, , drop = FALSE],
    cores = n.cores
  )
  # message("Finished mapping")
  # message(timestamp())
  peaks <- qtl2::find_peaks(out, gmap,
    drop = 1.5,
    threshold = thrA, thresholdX = thrX
  )
  peaks <- peaks %>%
    dplyr::select(-lodindex) %>%
    dplyr::rename(phenotype = lodcolumn, peak_chr = chr, peak_cM = pos)

  return(peaks)
}

batch_wrap <- function(tissue, exprZ_list, kinship_loco, qtlprobs,
                       covar_list, gmap, thrA, thrX, cores) {

  num.batches <- max(c(round(ncol(exprZ_list[[tissue]])/1000), 2))
  nn <- ncol(exprZ_list[[tissue]])
  ss <- round(seq(0, nn, length.out=num.batches + 1))

  # Calculate the maximum number of concurrent batches
  # Each batch should at least have 4 cores.
  max_concurrent_batches <- max(1, floor(cores / 4) )
  # if the #of batches > max_concurrent_batches adjust the cores
  if( num.batches > max_concurrent_batches){
    # Adjust the number of cores to use based on concurrency limit
    cores_to_use <- min(cores, max_concurrent_batches * 4)
    # get the cores to use per batch, minimum 4
    cores_per_batch <- max(4, floor(cores_to_use/max_concurrent_batches))
  } else{
    # can use all the cores
    cores_to_use <- cores
    # get the cores to use per batch
    cores_per_batch <- floor(cores_to_use/num.batches)
  }

  # Initialize an empty list to store the results
  all_results <- list()
  for (i in seq(1, num.batches, by = max_concurrent_batches)) {

    # Determine the current batch range
    current_batches <- i:min(i + max_concurrent_batches - 1, num.batches)

    # Run the current batch of tasks in parallel
    current_results <- foreach::foreach( i = current_batches)  %dopar% {
      peak_fun(i,
               ss,
               exprZ_list[[tissue]],
               kinship_loco[[tissue]],
               qtlprobs[[tissue]],
               covar_list[[tissue]],
               tissue = tissue,
               gmap = gmap,
               thrA = thrA,
               thrX = thrX,
               n.cores = cores_per_batch)
    }

    # Append the current results to the overall results
    all_results <- c(all_results, current_results)
  }


  peaks <- do.call("rbind", all_results)

  # message(paste0(colnames(peaks), sep = " "))

  tissue_peaks2 <- tibble::lst(tissue, peaks)

  return(tissue_peaks2)

}


get_cores <- function(){
  # get the total number of cores that are available for each OS
  if (Sys.info()['sysname'] == "Windows") {
    num_cores <- as.numeric(Sys.getenv("NUMBER_OF_PROCESSORS"))
    message(paste0("Working in Windows and there are ", num_cores, " cores in total."))
  }
  else if (Sys.info()['sysname'] == "Linux"){
    num_cores <- as.numeric(system("nproc", intern = TRUE))
    message(paste0("Working in Linux and there are ", num_cores, " cores in total."))
  }
  else if(Sys.info()['sysname'] == "Darwin" ){
    num_cores <- as.numeric(system("sysctl -n hw.ncpu", intern = TRUE))
    message(paste0("Working in MacOS and there are ", num_cores, " cores in total."))
  }
  else{
    num_cores <- 1
    message(paste0("Unknown OS using only 1 core."))
  }
  return(num_cores)
}

`%notin%` <- Negate(`%in%`)
