#' Plot an individual peak or Principal Component of a hotspot.
#'
#' @param mapping Mapping list from `mapQTL`.
#' @param feats List of phenotypes to include.
#' Individual/Candidate Mediator Only: Local to hotspot. PC: Targets of hotspot.
#' @param tbands List of transband locations (`bands.rna` list from
#'  `transbands` function).
#' @param chromosome Chromosome that the transband (hotspot) is present on.
#' @param tissue Tissue that the transband (hotspot) is present in.
#' @param candidateMed Dataframe of candidate mediator information
#'  (can be from annotations). Columns must include "id", "symbol",
#'   "start", "end".
#' @param hsNum If there are multiple hotspots on a chromosome,
#'  indicate which one. Default is 1.
#' @param color Plot color. Default is "#0073C2FF".
#' @param pc Logical. Plot the first principal component of the targets.
#' Default `FALSE`.
#' @param wag Logical. Include `tailWag` effects plot in final.
#' Default is `FALSE`.
#' @param sigLOD Significant LOD threshold to use for filtering phenotypes.
#' @param topPC Proportion of top hotspot contributers to be identified with
#'  Principal Component Analysis. Default is 0.5.
#' @param psave Logical. Save the plot as `.png`. Default `TRUE`.
#' @param pname File name to save plot (needs to end in `.png`).
#' @param outdir Directory to save plots. Default is `NULL`.
#' @param ... Additional arguments
#'
#' @return ggplot/patchwork object: Peak plot (with or without phenotype
#' effects). If plotting the Principal Component of the hotspot, also includes
#' a list of the top proportion of phenotypes contributing to the hotspot.
#' @export
#'
hsPeakPlot <- function(mapping, feats, tbands, chromosome, tissue,
                       candidateMed = NULL, hsNum = 1, color = "#0073C2FF",
                       pc = FALSE, wag = FALSE, sigLOD, topPC = 0.5,
                       psave = TRUE, pname = NULL, outdir = NULL, ...) {
  if (psave & is.null(outdir)) {
    stop("Plot to be saved, but no directory provided")
  }
  if (!pc & is.null(candidateMed)) {
    stop("If plotting a whole hotspot either a candidate mediator must be
         provided or ask to plot PC1")
  }
  if (psave & is.null(pname)) {
    if (!wag & pc) {
      pname <- paste0(tissue, "_chromosome_", chromosome, "_transband_", hsNum,
                      "_principal_component_peak.png")
    }
    if (wag & pc) {
      pname <- paste0(tissue, "_chromosome_", chromosome, "_transband_", hsNum,
                      "_principal_component_wag_peak.png")
    }
    if (!wag & !pc) {
      pname <- paste0(tissue, "_chromosome_", chromosome, "_transband_", hsNum,
                      "_peak.png")
    }
    if (wag & !pc) {
      pname <- paste0(tissue, "_chromosome_", chromosome, "_transband_", hsNum,
                      "_wag_peak.png")
    }
    message(paste0("Plot to be saved. Saving as ", pname, " in ", outdir))
  }

  ## Isolate the features to include
  hSpot <- tbands[[tissue]] |>
    dplyr::filter(chr == chromosome)

  hsStart <- hSpot$start[hsNum]
  hsEnd <- hSpot$end[hsNum]

  ## Calculate PC if including otherwise get rZ phenotypes
  if (pc) {
    counts <- mapping$expr_list[[tissue]][feats, , drop = FALSE]

    # return(counts)
    hsPCA <- stats::prcomp(x      = t(counts),
                           center = TRUE,
                           scale. = TRUE)

    pca_rz <- apply(hsPCA$x[,"PC1", drop=FALSE],2, rankZ)

    loadings_pc1 <- hsPCA$rotation[,"PC1", drop = FALSE] |> as.data.frame()
    cv <- sd(loadings_pc1$PC1) / mean(abs(loadings_pc1$PC1))

    if (cv < 0.05) {
      warning("Loadings for PC1 are very similar. Selecting top genes may not
              be meaningful.")
    }

    top_genes <- loadings_pc1 |>
      dplyr::slice_max(prop = topPC, order_by = abs(PC1)) |>
      rownames()
    # return(pca_rz)
  } else {
    feat_rz <- mapping$exprZ_list[[tissue]][,feats, drop = FALSE]
  }

  ## Perform phenotype scans
  if (pc) {
    scan_out <- qtl2::scan1(pheno     = pca_rz,
                            genoprobs = mapping$qtlprobs[[tissue]],
                            kinship   = mapping$kinship_loco[[tissue]],
                            addcovar  = mapping$covar_list[[tissue]])
    # return(scan_out)
  } else {
    scan_out <- qtl2::scan1(pheno     = feat_rz,
                            genoprobs = mapping$qtlprobs[[tissue]],
                            kinship   = mapping$kinship_loco[[tissue]],
                            addcovar  = mapping$covar_list[[tissue]])
  }

  if (!is.null(candidateMed)) {

    medID <- candidateMed$id
    medSym <- candidateMed$symbol
    medStart <- candidateMed$start
    medEnd <- candidateMed$end

    if (wag) {
      scan_og <- scan_out
    }

    if (!pc) {
      scan_out <- scan_out |>
        as.data.frame( ) |>
        tibble::rownames_to_column("marker") |>
        dplyr::rename(!!!purrr::set_names(medID, medSym)) |>
        dplyr::mutate( marker = dimnames(scan_out)[[1]]) |>
        dplyr::left_join(mapping$map_dat2)
    }
  }

  ## Perform effects scan if wanted
  if (wag) {
    ncores <- get_cores()
    cores_use <- min(8, ncores)
    if (pc & is.null(candidateMed)) {
      effects_out <- qtl2::scan1blup(
        pheno     = pca_rz,
        genoprobs = mapping$qtlprobs[[tissue]][, as.character(chromosome)],
        kinship   = mapping$kinship_loco[[tissue]][[as.character(chromosome)]],
        addcovar  = mapping$covar_list[[tissue]],
        cores     = cores_use)
    } else if (pc & !is.null(candidateMed)) {
      effects_out <- qtl2::scan1blup(
        pheno     = mapping$exprZ_list[[tissue]][, medID, drop = FALSE],
        genoprobs = mapping$qtlprobs[[tissue]][, as.character(chromosome)],
        kinship   = mapping$kinship_loco[[tissue]][[as.character(chromosome)]],
        addcovar  = mapping$covar_list[[tissue]],
        cores     = cores_use)
    } else {
      effects_out <- qtl2::scan1blup(
        pheno     = feat_rz[, medID],
        genoprobs = mapping$qtlprobs[[tissue]][, as.character(chromosome)],
        kinship   = mapping$kinship_loco[[tissue]][[as.character(chromosome)]],
        addcovar  = mapping$covar_list[[tissue]],
        cores     = cores_use)
    }
  }

  ## Set bounds for plotting
  chr_end <- max(mapping$map_dat2$pos_bp[which(mapping$map_dat2$chr ==
                                                 chromosome)])
  chr_start <- min(mapping$map_dat2$pos_bp[which(mapping$map_dat2$chr ==
                                                   chromosome)])
  hs_mid <- (hsStart + hsEnd) / 2

  minX <- max(chr_start, hs_mid - 25e6)
  maxX <- min(chr_end, hs_mid + 25e6)

  if ((maxX - minX) < 50e6) {
    if (minX == chr_start) {
      maxX <- min(chr_end, chr_start + 50e6)
    } else if (maxX == chr_end){
      minX <- max(chr_start, chr_end - 50e6)
    }
  }

  # message(minX)
  # message(maxX)

  ## Plotting
  if (pc & is.null(candidateMed)) {

    maxY <- ceiling(max(scan_out[, "PC1"] / 5)) * 5
    p1 <- scan_out |>
      tibble::as_tibble( rownames = "marker") |>
      dplyr::left_join(mapping$map_dat2) |>
      dplyr::filter( chr == chromosome) |>
      dplyr::mutate( type = paste0("Chr", chromosome, " PC1")) |>
      ggplot2::ggplot() +
      ggplot2::aes(
        x= pos_bp/1e06,
        y = PC1,
        col = type
      ) +
      ggplot2::geom_rect(  xmin = hsStart/1e06,
                           xmax = hsEnd/1e06,
                           ymin = 0,
                           ymax = maxY,
                           fill = "gray",
                           inherit.aes = FALSE,
                           alpha = 0.1,
                           show.legend = FALSE) +
      ggplot2::geom_line( size = 1.5) +
      ggplot2::theme_minimal(base_size = 18) +
      ggplot2::scale_color_manual(values = color) +
      ggplot2::xlab(paste0("Chr ",chromosome," location (Mbp)")) +
      ggplot2::ylab("LOD score") +
      ggplot2::ylim(0,maxY) +
      ggplot2::xlim(minX / 1e06, maxX / 1e06) +
      ggplot2::theme(legend.position = "none")
  } else if (pc & !is.null(candidateMed)) {
    maxY <- ceiling(max(scan_out[, "PC1"] / 5)) * 5
    p1 <- scan_out |>
      tibble::as_tibble( rownames = "marker") |>
      dplyr::left_join(mapping$map_dat2) |>
      dplyr::filter( chr == chromosome) |>
      dplyr::mutate(type = paste0("Chr", chromosome, " PC1")) |>
      ggplot2::ggplot() +
      ggplot2::aes(
        x= pos_bp/1e06,
        y = PC1,
        col = type
      ) +
      ggplot2::geom_rect(  xmin = hsStart/1e06,
                           xmax = hsEnd/1e06,
                           ymin = 0,
                           ymax = maxY,
                           fill = "gray",
                           inherit.aes = FALSE,
                           alpha = 0.1,
                           show.legend = FALSE) +
      ggplot2::geom_line( size = 1.5) +
      ggplot2::geom_segment(x = medStart/1e06,
                            xend = medEnd/1e06,
                            y = 0,
                            yend = 1,
                            size = 2,
                            ggplot2::aes(col = type)) +
      ggplot2::annotate("text", x= ((medStart/1e06) + (medEnd/1e06)) / 2,
                        y = -0.5, label = medSym, size = 4,
                        fontface = "italic") +
      ggplot2::theme_minimal(base_size = 18) +
      ggplot2::scale_color_manual(values = color) +
      ggplot2::xlab(paste0("Chr ",chromosome," location (Mbp)")) +
      ggplot2::ylab("LOD score") +
      ggplot2::ylim(-0.7, maxY) +
      ggplot2::xlim(minX / 1e06, maxX / 1e06) +
      ggplot2::theme(legend.position = "none")


  } else {
    maxY <- ceiling(max(scan_out[,c(feats[feats != medID], medSym)] / 5)) * 5

    # message(maxY)

    p1 <- scan_out |>
      dplyr::filter(chr == chromosome) |>
      tidyr::pivot_longer(cols = medSym,
                          names_to = "symbol",
                          values_to = "lod") |>
      ggplot2::ggplot() +
      ggplot2::aes(
        x = pos_bp/1e06,
        y = lod,
        col = symbol
      )+
      ggplot2::geom_rect(xmin = hsStart/1e06,
                         xmax = hsEnd/1e06,
                         ymin = 0,
                         ymax = maxY,
                         fill = "gray",
                         inherit.aes = FALSE,
                         alpha = 0.1,
                         show.legend = FALSE) +
      ggplot2::geom_line(size = 1.5, alpha = 0.8) +
      ggplot2::geom_segment(x = medStart/1e06,
                            xend = medEnd/1e06,
                            y = 0,
                            yend = 1,
                            size = 2,
                            ggplot2::aes(col = symbol)) +
      ggplot2::annotate("text", x= ((medStart/1e06) + (medEnd/1e06)) / 2,
                        y = -0.5, label = medSym, size =4,
                        fontface = "italic") +
      ggplot2::theme_minimal( base_size = 18) +
      ggplot2::scale_color_manual(values = color) +
      ggplot2::xlab(paste0("Chr ",chromosome ," location (Mbp)")) +
      ggplot2::ylab( "LOD score") +
      ggplot2::labs(col = "Gene") +
      ggplot2::ylim(-0.7, maxY) +
      ggplot2::xlim(minX / 1e06,
                    maxX / 1e06) +
      ggplot2::theme(legend.position = "none")
  }

  if (wag) {
    # return(scan_out)
    if (pc) {
      peaks <- qtl2::find_peaks(scan1_output = scan_out, map = mapping$pmap,
                                threshold = sigLOD, prob = 0.95)
    } else {
      peaks <- qtl2::find_peaks(scan1_output = scan_og, map = mapping$pmap,
                                threshold = sigLOD, prob = 0.95)
    }
    peaks <- peaks |>
      dplyr::select(-lodindex) |>
      dplyr::rename(phenotype = lodcolumn, peak_chr = chr, peak_bp = pos)

    if (!is.null(candidateMed) & !pc) {
      peaks <- peaks |>
        dplyr::filter(phenotype == medID)
    }

    correct_marker <- mapping$map_dat2 |>
      dplyr::filter(chr == chromosome) |>
      dplyr::slice_min(abs(pos_bp - as.numeric(peaks$peak_bp)), n = 1)

    peaks$marker <- correct_marker$marker
    effects_out <- effects_out[correct_marker$marker, , drop = FALSE]

    effects <- tibble::lst(peaks = setNames(list(peaks), tissue),
                           effects_blup = setNames(list(effects_out), tissue))

    if (!pc & !is.null(candidateMed)) {
      p2 <- tailWag(effects = effects, tissue = tissue, feat = medID,
                    chromosome = chromosome, psave = FALSE, symbol = FALSE, ...)
    } else {
      p2 <- tailWag(effects = effects, tissue = tissue, feat = "PC1",
                    chromosome = chromosome, psave = FALSE, symbol = FALSE,...)
    }
    pwork <- patchwork::wrap_plots(p1 + ggplot2::theme(legend.position = "top"),
                                   p2, nrow = 1, widths = c(1, 0.5))
  }

  if (psave) {
    if (wag) {
      ggplot2::ggsave(pname, plot = pwork, device = "png", path = outdir,
                      bg = "white", width = 7, height = 5, units = "in")
    } else {
      ggplot2::ggsave(pname, plot = p1, device = "png", path = outdir,
                      bg = "white", wdth = 5, height = 5, units = "in")
    }
  }
  if (wag) {
    if(pc) {
      return(tibble::lst(pca_rz, top_genes, pwork, peaks))
    } else {
      return(pwork)
    }
  } else if (pc) return(tibble::lst(pca_rz, top_genes, p1))
  else return(p1)
}
