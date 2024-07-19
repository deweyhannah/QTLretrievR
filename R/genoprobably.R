## usethis namespace: start
#' @export
#'
#' @importFrom abind abind
## usethis namespace: end
## Previously called tsv2genoprobs
genoprobably <- function(outfile = "./gbrs_interpolated_genoprobs.RData", gbrsFileLoc, tissues = c(), gridfile = "/projects/compsci/omics_share/mouse/GRCm39/supporting_files/emase_gbrs/rel_2112_v8/ref.genome_grid.GRCm39.tsv"){

  ## Check that interpolated genoprobs doesn't already exist
  gbrs.interp <- outfile
  if (file.exists(gbrs.interp)) {
    stop(paste0(outfile, " allready exists. Please delete or use new file name."))
  }

  ## Check that there is at least one sample/tissue
  filenames <- c(list.files(gbrsFileLoc,
                            pattern='.*.gbrs.interpolated.genoprobs.tsv$', full.names = TRUE, recursive = T) )
  if(length(filenames) == 0){
    stop(paste0(gbrsFileLoc, " does not contain any interpolated genoprobs"))
  }
  if(length(filenames)/length(tissues) < 1 & length(tissues)!=0){
    stop(paste0(gbrsFileLoc, " does not contian enough files for the tissues indicated"))
  }

  ## Separate files by tissue type - tissue names provided should match how they are named in the samples
  tissue_files <- list()
  if(length(tissues) > 0){
    for(tissue in tissues){
      tissue_files[[tissue]] <- filenames[grepl(tissue, filenames)]
    }
  }

  if(length(tissues) == 0){
    tissue_files[["a"]] <- filenames
  }

  ## Split file names into sample names for each tissue type
  tissue_samps <- list()
  for(tissue in tissues){
    samples <- unlist(strsplit(basename(tissue_files[[tissue]]), ".", fixed = TRUE))[[1]]
    for(i in 2:length(tissue_files[[tissue]])){
      samples[i] <- unlist(strsplit(basename(tissue_files[[tissue]])[i], ".", fixed = TRUE))[[1]]
    }
    tissue_samps[[tissue]] <- samples
  }

  ## Load in genotype tsvs for each sample/tissue
  genoprobs_list <- list()
  for(tissue in tissues){
    genoprobs_list[[tissue]] <- list()
    for (i in 1:length(tissue_samps[[tissue]])) {
      tsvfile <- tissue_files[[tissue]][i]
      samp <- tissue_samps[[tissue]][i]
      message(i)
      genoprobs_list[[tissue]][[samp]] <- readr::read_tsv(tsvfile, col_names = LETTERS[1:8], show_col_types = FALSE, skip = 1)
    }
  }

  ## Organize genoprobs to 3d structure
  tsv_probs <- list()
  for(tissue in tissues){
    tsv_probs[[tissue]] <- abind::abind(genoprobs_list[[tissue]], along = 3)
    tsv_probs[[tissue]] <- aperm(tsv_probs[[tissue]], c(3, 2, 1))
  }

  ## Attach marker names to probs
  gridmap <- readr::read_tsv(gridfile, col_types='ccddd')
  for(tissue in tissues){
    dimnames(tsv_probs[[tissue]])[[3]] <- gridmap$marker
  }

  ## Rename
  probs_list <- tsv_probs

  ## Save to outfile
  saveRDS(probs_list, gbrs.interp)

  return(probs_list)
}
