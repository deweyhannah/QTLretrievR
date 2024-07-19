## usethis namespace: start
#' @export
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @import tibble
#' @import tidyr
## usethis namespace: end
mapQTL <- function(outdir, peaks_out, map_out, genoprobs, samp_meta, expr_mats, covar_factors, n.cores = 4, gridfile = "/projects/compsci/omics_share/mouse/GRCm39/supporting_files/emase_gbrs/rel_2112_v8/ref.genome_grid.GRCm39.tsv", localRange = 10e6,
                    biomart, samp_excl = c()){
  ## Expression Matrices should be listed in the same order as tissues were for tsv2genoprobs call
  ## Load probs
  if(is.list(genoprobs)){
    tmp_probs <- check_data(genoprobs)
  }
  if(is.character(genoprobs)){
    tmp_probs <- check_data(paste0(outdir,"/",genoprobs), "genoprobs")
  }
  probs_list <- tmp_probs
  rm(tmp_probs)

  ## Check inputs

  if(length(expr_mats) == 0){
    stop("Please provide at least one expression matrix")
  }
  if(length(names(probs_list)) > length(expr_mats)){
    stop("More probabilites found than expression matrices. Are you missing any expression inputs?")
  }
  if(length(names(probs_list)) < length(expr_mats)){
    stop("More expression matrices than probabilites found. Are you using the correct probabilities?")
  }
  if(length(covar_factors) == 0){
    stop("Please provide at least one covariate for model")
  }

  ## Modify Probs and Determine Kinship
  qtlprobs <- list()
  kinship_loco <- list()
  for(tissue in names(probs_list)){
    message(tissue)
    qtlprobs[[tissue]] <- probs_3d_to_qtl2(probs_list[[tissue]])
    kinship_loco[[tissue]] <- qtl2::calc_kinship(qtlprobs[[tissue]], "loco", cores = n.cores)
  }

  ## Create maps
  grid_map <- read.delim(gridfile, stringsAsFactors = F, row.names = 1)
  map_dat <- grid_map[,c("chr", "cM")]

  markers <- dimnames(probs_list[[1]])[[3]]
  map_dat <-map_dat[markers, ]
  # assert_that(are_equal(rownames(map_dat), markers))
  gmap <- split_map(map_dat)
  map_dat$marker <- rownames(map_dat)

  message("gmap complete")

  map_dat2 <- map_dat %>%
    tidyr::separate(marker, into=c('chrom', 'pos_bp'), convert=T, remove=F) %>%
    dplyr::mutate(n=1:dplyr::n()) %>%
    tibble::as_tibble()
  pmap <- split_map(dplyr::select(map_dat2, marker,
                                  chr, pos_bp) %>%
                      as.data.frame() %>%
                      tibble::remove_rownames() %>%
                      tibble::column_to_rownames('marker'))

  message("pmap complete")

  ## Load expression matrices
  expr_list <- list()
  for(i in 1:length(expr_mats)){
    expr <- expr_mats[i]
    tissue <- names(probs_list)[[i]]
    expr_list[[tissue]] <- as.matrix(read.delim(expr, row.names = 1, header = T, stringsAsFactors = F))
  }

  message("expression loaded")

  ## Sample Details
  sample_details <- read.delim(samp_meta, stringsAsFactors = F, header = T)
  if(!all(covar_factors %in% colnames(sample_details))){
    stop("Chosen factors are not in sample metadata. Please check factors and sample metadata for missing or misspelled elements")
  }
  for(fact in covar_factors){
    sample_details[,fact] <- as.factor(sample_details[,fact])
  }

  ## Reorganize and calculate rankZ for expression matrices
  exprZ_list <- list()
  for(tissue in names(expr_list)){
    message(tissue)
    samps_keep <- rownames(probs_list[[tissue]])[which(rownames(probs_list[[tissue]]) %notin% samp_excl)]
    message(length(samp_excl))
    message(length(samps_keep))
    message(length(colnames(expr_list[[tissue]])))
    message(length(intersect(colnames(expr_list[[tissue]]), samps_keep)))
    expr_list[[tissue]] <- expr_list[[tissue]][, samps_keep, drop = FALSE]
    exprZ_list[[tissue]] <- apply(expr_list[[tissue]], 1, rankZ)
  }

  std <- standardize(expr_list, exprZ_list, sample_details, tissues = names(expr_list))
  expr_list <- std$expr
  exprZ_list <- std$exprZ
  tissue_samp <- std$tissue_samp

  message("rankZ normalized")

  ## Calculate covariate matrices
  covar_list <- list()
  for(tissue in names(tissue_samp)){
    covar_list[[tissue]] <- model.matrix( formula(paste0("~",paste0(covar_factors, collapse = "+"))), data = tissue_samp[[tissue]])
    rownames(covar_list[[tissue]]) <- tissue_samp[[tissue]]$ID
    covar_list[[tissue]] <- covar_list[[tissue]][, colnames(covar_list[[tissue]]) != "(Intercept)", drop = FALSE]
    for(fact in covar_factors){
      if(length(levels(tissue_samp[[tissue]][,fact])) == 2){
        colnames(covar_list[[tissue]])[grepl(fact, colnames(covar_list[[tissue]]))] <- fact
      }
    }
  }

  message("covariates calculated")

  outfile <- paste0(outdir, "/", map_out)

  maps_list <- list(qtlprobs, covar_list, expr_list, exprZ_list, kinship_loco, gmap, map_dat2,pmap, tissue_samp)
  names(maps_list) <- c("qtlprobs", "covar_list", "expr_list", "exprZ_list", "kinship_loco", "gmap", "map_dat2", "pmap", "tissue_samp")
  saveRDS(maps_list, file = outfile)

  ## Run Batchmap

  message("calculating peaks")

  ## Add functionality to save out files if wanted -- probably in the background functions file.
  peaks_list <- list()
  for(tissue in names(exprZ_list)){
    num.batches <- max(c(round(ncol(exprZ_list[[tissue]])/1000), 2))
    message("Mapping ", ncol(exprZ_list[[tissue]]), " ", tissue, " gene expression levels. Running in ", num.batches, " batches")
    peaks_list[[tissue]] <- batchmap(num.batches, exprZ_list[[tissue]], kinship_loco[[tissue]], qtlprobs[[tissue]], covar_list[[tissue]], tissue, gmap = gmap, n.cores = n.cores)
  }

  peaks_list <- annotatePeaks(maps_list, peaks_list, biomart, localRange)

  outfile <- paste0(outdir,"/", peaks_out)
  saveRDS(peaks_list, file=outfile)

  map_peaks <- list(maps_list, peaks_list)
  names(map_peaks) <- c("maps_list", "peaks_list")
  return(map_peaks)

}
