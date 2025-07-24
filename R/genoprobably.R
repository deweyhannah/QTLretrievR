#' Convert GBRS tsv genome probabilities to genoprobs format.
#' @description
#' Converts GBRS-formatted TSV genotype probability files into a 3D array format compatible with `qtl2`, organized by tissue and sample, optionally saves the result as an RDS file
#'
#'
#' @param outfile File path to save the resulting genopbrobs list. Defaults to "gbrs_interpolated_genoprobs.rds".
#' @param gbrsFileLoc File path to where the GBRS files are located.
#' @param tissues Character vector of tissue names. These should match substrings in the GBRS file names. if empty, defaults fo "a".
#' @param gridFile Either a data frame or file path to the genome grid used to assign marker names. Default is 75K grid loaded with the package.
#' @param save Character. Determines output behavior: "sr" to save and return, "so" to save only, "ro" to return only. Default is "sr".
#'
#' @return A named list of 3D genotype probability arrays (one per tissue), formatted for use with `qtl2`.
#' @export
#'
#' @importFrom abind abind
#'
genoprobably <- function(outfile = "./gbrs_interpolated_genoprobs.rds",
                         gbrsFileLoc,
                         tissues = c(),
                         gridFile = gridfile,
                         save = "sr") {
  ## Check that interpolated genoprobs doesn't already exist
  gbrs.interp <- outfile
  if (file.exists(gbrs.interp)) {
    stop(paste0(outfile, " allready exists. Please delete or use new file name."))
  }

  ## Locate all interpolated genoprob TSV files in the specified directory
  filenames <- c(list.files(gbrsFileLoc,
    pattern = ".*.gbrs.interpolated.genoprobs.tsv$", full.names = TRUE, recursive = T
  ))
  if (length(filenames) == 0) {
    stop(paste0(gbrsFileLoc, " does not contain any interpolated genoprobs"))
  }
  if (length(tissues) == 0) {
    tissues <- "a"
  }
  ## Ensure there are enough files to match the number of tissues provided
  if (length(filenames) / length(tissues) < 1 & length(tissues) != 0) {
    stop(paste0(gbrsFileLoc, " does not contian enough files for the tissues indicated"))
  }

  ## Separate files by tissue type - tissue names provided should match how they are named in the samples
  tissue_files <- list()
  if (length(tissues) > 0) {
    for (tissue in tissues) {
      tissue_files[[tissue]] <- filenames[grepl(tissue, filenames)]
    }
  }

  ## Extract sample names from file names for each tissue
  tissue_samps <- list()
  for (tissue in tissues) {
    samples <- unlist(strsplit(basename(tissue_files[[tissue]]), ".", fixed = TRUE))[[1]]
    for (i in 2:length(tissue_files[[tissue]])) {
      samples[i] <- unlist(strsplit(basename(tissue_files[[tissue]])[i], ".", fixed = TRUE))[[1]]
    }
    tissue_samps[[tissue]] <- samples
  }

  ## Read each TSV file into a list of matrices, organized by tissue and sample
  genoprobs_list <- list()
  for (tissue in tissues) {
    genoprobs_list[[tissue]] <- list()
    for (i in 1:length(tissue_samps[[tissue]])) {
      tsvfile <- tissue_files[[tissue]][i]
      samp <- tissue_samps[[tissue]][i]
      # message(i)
      genoprobs_list[[tissue]][[samp]] <- readr::read_tsv(tsvfile, col_names = LETTERS[1:8], show_col_types = FALSE, skip = 1)
    }
  }

  ## Combine individual sample matrices into a 3D array and reorder dimensions
  tsv_probs <- list()
  for (tissue in tissues) {
    tsv_probs[[tissue]] <- abind::abind(genoprobs_list[[tissue]], along = 3)
    tsv_probs[[tissue]] <- aperm(tsv_probs[[tissue]], c(3, 2, 1))
  }

  ## Assign marker names to the third dimension of the 3D array using the genome grid
  if( is.data.frame(gridFile)) {
    gridmap <- gridFile
    for (tissue in tissues) {
      dimnames(tsv_probs[[tissue]])[[3]] <- rownames(gridmap)
    }
  }
  if( is.character(gridFile)) {
    gridmap <- readr::read_tsv(gridFile, col_types = "ccddd")
    for (tissue in tissues) {
      dimnames(tsv_probs[[tissue]])[[3]] <- gridmap$marker
    }
  }

  ## Rename
  # probs_list <- tsv_probs
  probs_list <- list()
  for (tissue in names(tsv_probs)) {
    message(tissue)
    probs_list[[tissue]] <- probs_3d_to_qtl2(tsv_probs[[tissue]])
  }

  ## Save to outfile
  if(save %in% c("sr", "so")) {
    saveRDS(probs_list, gbrs.interp)
  }

  if(save %in% c("sr","ro")) {
    return(probs_list)
  }
}
