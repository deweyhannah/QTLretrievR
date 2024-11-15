#' Convert and process MUGA probabilities to qtl2 format
#'
#' @param type string indicating GigaMUGA (GM) or MegaMUGA (MM), default is "GM"
#' @param covarLoc location of covariate file
#' @param covar_file covariate file including at minimum the sex (sex) and generation (ngen) of each sample, this needs to be a .csv file
#' @param i.files either a string of the directory where the chromosome specific genotype files are or a list of final report files to process - if passing the final report files they need to be either unzipped or in .gz format
#' @param genoPrefix prefix for the chromosome specific genotype files (excluding "_geno")
#' @param probsOut file name to save probabilities, default is "muga_interpolated_genoprobs.rds". probs will save to directory that mugaprobs is called from.
#'
#' @return none
#' @export
#'
#' @importFrom httr GET
#' @importFrom httr write_disk
#' @importFrom data.table fread
#' @importFrom qtl2 calc_genoprob
#'
mugaprobs <- function(type = "GM", covarLoc, covar_file, i.files, genoPrefix = "gm4qtl2", probsOut = "muga_interpolated_geonprobs.rds") {
  ## Confirm type and set url to pull info from
  if(type == "GM") {
    message("Using GigaMUGA markers for calculating probabilities")
    url <- "https://figshare.com/ndownloader/files/40233652"
    code_file <- "GM/GM_allelecodes.csv"
  }
  if(type == "MM") {
    message("Using MegaMUGA markers for calculating probabilities")
    url <- "https://figshare.com/ndownloader/files/9311164"
    code_file <- "MM/MM_allelecodes.csv"
  }

  ## Pull relevant files for MUGA
  temp_dir <- tempdir()
  temp_zip <- tempfile(fileext = ".zip")

  response <- httr::GET(url, httr::write_disk(temp_zip, overwrite = TRUE))
  unzip(temp_zip, exdir = temp_dir, overwrite = T)

  file.copy(paste0(covarLoc, "/" , covar_file), to = temp_dir)

  ## If passed files are not a directory process final reports
  if(!dir.exists(i.files)) {
    message("processing final reports")
    process_reports(paste0(temp_dir,"/",code_file), i.files, genoPrefix, temp_dir)
  }

  ## Get current directory to go back after processing and move to temp directory
  ogDir <- getwd()
  setwd(temp_dir)

  ## Bring chromosome specific genotype files to temp directory
  if(dir.exists(i.files)) {
    files_to_copy <- list.files(i.files, pattern = paste0(genoPrefix,"_geno"), full.names = T)
    file.copy(files_to_copy, to = temp_dir)
  }

  ## Write control file
  message("writing control file")
  chr <- c(1:19,"X")
  write_control_file(paste0("forqtl2.json"),
                     crosstype = "do",
                     description = "DO Project",
                     founder_geno_file = paste0(type, "/", type, "_foundergeno", chr, ".csv"),
                     founder_geno_transposed = TRUE,
                     gmap_file = paste0(type, "/", type, "_gmap", chr, ".csv"),
                     pmap_file = paste0(type, "/", type, "_pmap", chr, ".csv"),
                     geno_file = paste0(genoPrefix, "_geno", chr, ".csv"),
                     geno_transposed = TRUE,
                     geno_codes = list(A=1, H=2, B=3),
                     xchr = "X",
                     covar_file = covar_file,
                     sex_covar = "sex",
                     sex_codes = list(F="Female", M="Male", f = "Female", m = "Male"),
                     crossinfo_covar = "ngen",
                     overwrite = TRUE)

  ctrl <- read_cross2("forqtl2.json")
  map <- ctrl$pmap
  message("calculating genoprobs")
  pr <- qtl2::calc_genoprob(ctrl, map, error_prob = 0.002, cores = 4)

  saveRDS(pr, paste0(ogDir,"/", probsOut))
  setwd(ogDir)
}


#' @import qtl2convert
#' @importFrom data.table fread
process_reports <- function(codefile, ifiles, ostem, dirOut) {
  codes <- read.csv(codefile, comment.char = "#")

  full_geno <- NULL
  cXint <- cYint <- NULL

  for(ifile in ifiles) {
    rezip <- FALSE
    if(!file.exists(ifile)) {
      system(paste("gunzip", ifile))
      rezip <- TRUE
    }

    message(paste0("reading: ", ifile))
    # g <- read_reports(ifile)
    g <- data.table::fread(ifile, skip = 8, data.table = F)
    g <- g[g[,"SNP Name"] %in% codes[,"marker"],]

    # NOTE: may need to revise the IDs in the 2nd column
    samples <- unique(g[,"Sample ID"])

    # matrix to contain the genotypes
    geno <- matrix(nrow=nrow(codes), ncol=length(samples))
    dimnames(geno) <- list(codes[,"marker"], samples)

    # fill in matrix
    for(i in seq(along=samples)) {
      if(i==round(i,-1)) cat(" --Sample", i, "of", length(samples), "\n")
      wh <- (g[,"Sample ID"]==samples[i])
      geno[g[wh,"SNP Name"],i] <- paste0(g[wh,"Allele1 - Forward"], g[wh,"Allele2 - Forward"])
    }

    geno <- qtl2convert::encode_geno(geno, as.matrix(codes[,c("A","B")]))

    if(is.null(full_geno)) {
      full_geno <- geno
    } else {
      # if any columns in both, use those from second set
      full_geno <- qtl2convert::cbind_smother(full_geno, geno)
    }

    if(rezip) {
      system(paste("gzip", ifile))
    }
  }

  # write data to chromosome-specific files
  message("writing chromosome-specific files")
  for(chr in c(1:19,"X","Y","M")) {
    mar <- codes[codes$chr==chr,"marker"]
    g <- full_geno[mar,]
    qtl2convert::write2csv(cbind(marker=rownames(g), g),
                           paste0(dirOut, "/", ostem, "_geno", chr, ".csv"),
                           paste0(ostem, " genotypes for chr ", chr),
                           overwrite=TRUE)
  }
}

# read_reports <- function(file_path) {
#   lines <- readLines(file_path)
#   # Find the first non-comment line
#   data_start <- which(grepl("^SNP", lines))[1]
#   # Read the data starting from the first non-comment line
#   # # data <- read.delim(text = paste(lines[data_start:length(lines)], collapse = "\n"), header = TRUE, check.names = F)
#   # data <- readr::read_delim(file_path, delim = "\t", skip = data_start - 1, col_names = TRUE, comment = "#", show_col_types = FALSE)
#
#   return(data)
# }
