#' Background functions for QTLretrievR that will be called by the primary functions
#' @importFrom dplyr rename select arrange
#' @importFrom DT datatable
#' @importFrom BiocParallel SerialParam MulticoreParam bpnworkers
#' @importFrom stats setNames
#' @import qtl2
#' @importFrom intermediate mediation.scan
#' @importFrom stats qnorm
#'

#' @param map grid map of markers, chromosomes, and positions of markers
#' @param chr_names vector of chromosome names
split_map <- function(map, chr_names = NULL) {
  map <- reorder_map_table(map, chr_names = chr_names)
  pos <- as.numeric(map[, 2])
  chr <- map[, 1]
  uchr <- unique(chr)
  names(pos) <- rownames(map)
  lapply(split(pos, factor(chr, uchr)), sort)
}

reorder_map_table <- function(map_tab, chr_col = 1,
                              pos_col = 2, chr_names = NULL) {
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
  chroms <- vapply(strsplit(markers, "_"), function(x) x[[1]], character(1))
  newprobs <- vector("list", length(uchroms))
  names(newprobs) <- uchroms
  for (chrom in uchroms) newprobs[[chrom]] <- probs[, , chroms == chrom]
  attr(newprobs, "is_x_chr") <- c(rep(FALSE, length(uchroms) - 1), TRUE)
  attr(newprobs, "crosstype") <- "do"
  attr(newprobs, "alleles") <- c("A", "B", "C", "D", "E", "F", "G", "H")
  attr(newprobs, "alleleprobs") <- TRUE
  class(newprobs) <- c("calc_genoprob", "list")
  newprobs
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

standardize <- function(expr, details, tissues = c()) {
  tissue_samp <- list()
  for (tissue in tissues) {
    tissue_samp[[tissue]] <- details[which(details$ID %in%
                                             colnames(expr[[tissue]])) ,]
  }
  std <- tibble::lst(expr, tissue_samp)
  # names(std) <- c("expr","exprZ","tissue_samp")
  return(std)
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

peak_fun <- function(i,ss, exprZ, kinship_loco, genoprobs, covar, xcovar, tissue, gmap,
                     thrA = 5, thrX = 5, n.cores = 4, phys) {

  # message("Started mapping")
  # message(timestamp())
  start <- ss[i] + 1
  end <- ss[i + 1]
  # out <- qtl2::scan1(
  #   genoprobs,
  #   exprZ[,start:end,drop=F],
  #   kinship_loco,
  #   addcovar = covar[, , drop = FALSE],
  #   cores = n.cores
  # )

  # return(tibble::lst(genoprobs, exprZ, kinship_loco, covar))

  # message(sprintf("Running batch %d: cols %d-%d with %d cores", i, start, end, n.cores))
  # message(timestamp())

  out <- tryCatch({
    qtl2::scan1(
      genoprobs,
      exprZ[, start:end, drop = FALSE],
      kinship_loco,
      addcovar = covar[, , drop = FALSE],
      Xcovar = xcovar[, , drop = FALSE],
      cores = n.cores
    )
  }, error = function(e) {
    message("ERROR in batch: ", i)
    message("Slice: ", start, "-", end)
    message("Expr dims: ", paste(dim(exprZ[, start:end, drop=FALSE]), collapse=" x "))
    message("Any NA in expr: ", any(is.na(exprZ[, start:end])))
    message("Covar dims: ", paste(dim(covar), collapse=" x "))
    message("Xcovar dims: ", if (!is.null(xcovar)) paste(dim(xcovar), collapse=" x ") else "NULL")
    stop(e)
  })

  # message("Finished mapping")
  # message(timestamp())

  peaks <- qtl2::find_peaks(scan1_output = out,
                            map          = gmap,
                            prob         = 0.95,
                            threshold    = thrA,
                            thresholdX   = thrX
  )
  if (phys) {
    peaks <- peaks |>
      dplyr::select(-lodindex) |>
      dplyr::rename(phenotype = lodcolumn, peak_chr = chr, peak_bp = pos)
  } else {
    peaks <- peaks |>
      dplyr::select(-lodindex) |>
      dplyr::rename(phenotype = lodcolumn, peak_chr = chr, peak_cM = pos)
  }

  return(peaks)
}

batch_wrap <- function(tissue, exprZ_list, kinship_loco, qtlprobs,
                       covar_list, xcovar_list, gmap, thrA, thrX, BPPARAM, phys, min_cores = 4) {

  workers <- BiocParallel::bpnworkers(BPPARAM)

  num.batches <- max(c(round(ncol(exprZ_list[[tissue]])/1000), 2))
  nn <- ncol(exprZ_list[[tissue]])
  ss <- round(seq(0, nn, length.out=num.batches + 1))

  # Calculate the maximum number of concurrent batches
  # Each batch should at least have 4 cores.
  max_concurrent_batches <- max(1, floor(workers / min_cores) )
  # if the #of batches > max_concurrent_batches adjust the cores
  if( num.batches > max_concurrent_batches){
    # Adjust the number of cores to use based on concurrency limit
    cores_to_use <- min(workers, max_concurrent_batches * min_cores)
    # get the cores to use per batch, minimum 4
    cores_per_batch <- max(min_cores,
                           floor(cores_to_use/max_concurrent_batches))
  } else{
    # can use all the cores
    cores_to_use <- workers
    # get the cores to use per batch
    cores_per_batch <- floor(cores_to_use/num.batches)
  }

  # Initialize an empty list to store the results
  all_results <- list()
  for (i in seq(1, num.batches, by = max_concurrent_batches)) {

    # Determine the current batch range
    current_batches <- i:min(i + max_concurrent_batches - 1, num.batches)

    # Run the current batch of tasks in parallel

    current_results <- BiocParallel::bplapply(current_batches,
                                              FUN = peak_fun,
                                              ss = ss,
                                              exprZ = exprZ_list[[tissue]],
                                              kinship_loco = kinship_loco[[tissue]],
                                              genoprobs = qtlprobs[[tissue]],
                                              covar = covar_list[[tissue]],
                                              xcovar = xcovar_list[[tissue]],
                                              tissue = tissue,
                                              gmap = gmap,
                                              thrA = thrA,
                                              thrX = thrX,
                                              n.cores = cores_per_batch,
                                              phys = phys,
                                              BPPARAM = BPPARAM)

    # Append the current results to the overall results
    all_results <- c(all_results, current_results)
  }



  peaks <- do.call("rbind", all_results)

  tissue_peaks2 <- tibble::lst(tissue, peaks)

  return(tissue_peaks2)
}

## Changes to this function made by Sam Widmayer (sam-widmayer)
get_cores <- function(){
  ntasks        <- Sys.getenv("SLURM_NTASKS")
  cpus_per_task <- Sys.getenv("SLURM_CPUS_PER_TASK")
  job_cpus      <- Sys.getenv("SLURM_JOB_CPUS_PER_NODE")

  if (ntasks != "" && cpus_per_task != "") {
    num_cores <- as.numeric(ntasks) * as.numeric(cpus_per_task)
    message(paste0("SLURM job detected (ntasks * cpus-per-task), using ",
                   num_cores, " cores."))
  } else if (cpus_per_task != "") {
    # Most common SLURM case: --cpus-per-task without explicit --ntasks
    num_cores <- as.numeric(cpus_per_task)
    message(paste0("SLURM job detected (cpus-per-task), using ",
                   num_cores, " cores."))
  } else if (job_cpus != "") {
    # SLURM_JOB_CPUS_PER_NODE can be formatted as "16(x2)" — take the number
    num_cores <- as.numeric(gsub("\\(.*", "", job_cpus))
    message(paste0("SLURM job detected (job-cpus-per-node), using ",
                   num_cores, " cores."))
  } else if (Sys.info()['sysname'] == "Windows") {
    num_cores <- as.numeric(Sys.getenv("NUMBER_OF_PROCESSORS"))
    message(paste0("Working in Windows and there are ", num_cores,
                   " cores in total."))
  } else if (Sys.info()['sysname'] == "Linux") {
    num_cores <- as.numeric(system("nproc", intern = TRUE))
    message(paste0("Working in Linux and there are ", num_cores,
                   " cores in total."))
  } else if (Sys.info()['sysname'] == "Darwin") {
    num_cores <- as.numeric(system("sysctl -n hw.ncpu", intern = TRUE))
    message(paste0("Working in MacOS and there are ", num_cores,
                   " cores in total."))
  } else {
    num_cores <- 1
    message("Unknown OS, using only 1 core.")
  }
  return(num_cores)
}

`%notin%` <- Negate(`%in%`)

filter_peaks <- function(peaks_df, bands_df) {
  purrr::map_dfr(seq_len(nrow(bands_df)), function(i) {
    band <- bands_df[i, ]
    peaks_df |>
      dplyr::filter(
        as.character(peak_chr) == as.character(band$chr),
        dplyr::between(peak_bp, band$start, band$end)
      )
  })
}

## This is essentially the same as qtl2convert::cbind_smother with one
## difference to fix an issue with replacement

cbind_smother_fix <- function(mat1, mat2)
  {
    cn1 <- colnames(mat1)
    cn2 <- colnames(mat2)
    if(is.null(cn1) || is.null(cn2)) {
      stop("Need colnames for both mat1 and mat2")
    }

    m_col <- (cn2 %in% cn1)

    if(any(m_col)) {
      rn1 <- rownames(mat1)
      rn2 <- rownames(mat2)
      if(is.null(rn1) || is.null(rn2)) {
        stop("Need rownames for both mat1 and mat2")
      }
      m_row <- match(rn1[rn1 %in% rn2], rn2)

      mat1 <- qtl2::cbind_expand(mat1[,!(cn1 %in% cn2[m_col]),drop=FALSE], mat2)
    }
    else {
      mat1 <- qtl2::cbind_expand(mat1, mat2)
    }

    mat1
  }
