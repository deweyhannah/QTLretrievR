#' Plot peaks associated with specific genes.
#'
#' @param mapping Mapping list from `mapQTL`, or full path to `.rds` containing
#' one.
#' @param tissue Tissue to derive plots from
#' @param pheno Phenotype to plot
#' @param pop One of `c("do", "cc", "other")` to indicate founder population.
#' Default is "`do`".
#' @param psave Logical. Save the plot as `.png`. Default `TRUE`.
#' @param pname File name to save plot (needs to end in `.png`).
#' @param outdir Directory to save plots. Default is `NULL`.
#' @param chromosome Chromosome that the peaks is present on.
#' @param effects Logical. Plotting effects on chromosome. Default is `FALSE`.
#' @param founders If `pop == "other"`, list of founders in haplotype order.
#' @param palette Founder color map (> 8 founders).
#'
#' @return Plot of phenotype specific peaks, or phenotype specific peak with
#'  effects
#' @export
#'
#' @importFrom qtl2 scan1 scan1blup plot_coefCC plot_coef
#' @importFrom grDevices png dev.off
#' @importFrom graphics plot.new legend
#'
peak_plot <- function(mapping, tissue, pheno, pop = "do", outdir = NULL,
                      pname = NULL, psave = TRUE, chromosome = NULL,
                      effects = FALSE, founders = NULL, palette = NULL) {
  if (is.list(mapping)) {
    tmp_map <- check_data(mapping)

  }
  if (is.character(mapping)) {
    tmp_map <- check_data(mapping)
  }
  if (!is.null(founders) & length(founders) > 8 & is.null(palette)) {
    stop(paste0(length(founders), " founders detected, please provide a palette
                that contains at least that many colors"))
  }

  if (psave & is.null(outdir)) {
    stop("Plot to be saved, but no directory provided")
  }
  if (psave & !is.null(pname)) {
    message(paste0("Plot to be saved. Saving as ", pname, " in ", outdir))
  }
  if (psave & is.null(pname)) {
    if (effects) {
      temp_name <- paste0(tissue,"_",pheno,"_effects.png")
    } else {
      temp_name <- paste0(tissue,"_",pheno,"_LOD.png")
    }
    message(paste0("Plot to be saved, but name not provided. Saving as ",
                   temp_name, " in ", outdir))
  }

  if (!is.null(tmp_map)) {
    exprZ_list <- tmp_map$exprZ_list
    covar_list <- tmp_map$covar_list
    expr_list <- tmp_map$expr_list
    gmap <- tmp_map$gmap
    kinship_loco <- tmp_map$kinship_loco
    map_dat2 <- tmp_map$map_dat2
    pmap <- tmp_map$pmap
    qtlprobs <- tmp_map$qtlprobs
    tissue_samp <- tmp_map$tissue_samp
  }
  rm(tmp_map)
  rm(mapping)

  if (is.null(chromosome) & effects) {
    stop("Peak chromosome needed to produce effects plot")
  }

  avail_cores <- get_cores()

  scan1_out <- qtl2::scan1(genoprobs = qtlprobs[[tissue]],
                           pheno     = exprZ_list[[tissue]][, pheno,
                                                            drop = FALSE],
                           kinship   = kinship_loco[[tissue]],
                           addcovar  = covar_list[[tissue]],
                           cores     = min(4, avail_cores))

  if (effects) {
    message("calculating effects")
    c2eff <- qtl2::scan1blup(
      genoprobs = qtlprobs[[tissue]][,as.character(chromosome)],
      pheno     = exprZ_list[[tissue]][, pheno, drop = FALSE])


    if (pop %in% c("do", "cc")) {
      p <- qtl2::plot_coefCC(x            = c2eff,
                             map          = pmap,
                             scan1_output = scan1_out,
                             bgcolor      = "gray95")
      if (psave) {
        if (is.null(pname)) {
          pname <- paste0(tissue, "_", pheno, "_effects.png")
        }
        png(paste0(outdir, "/", pname), width = 1024, height = 1024)
        print(p)
        dev.off()

        legend_plot <- paste0(outdir,"/CC_effects_legend.png")
        if (!file.exists(legend_plot)) {
          png(legend_plot, bg = "transparent")
          plot.new()
          legend("center", c("AJ","B6","129","NOD","NZO","CAST","PWK","WSB"),
                 lty = 1, col = qtl2::CCcolors, bty = "n", ncol = 4, lwd = 1.3)
          dev.off()
        }

      }
    }

    if (pop %notin% c("do", "cc")) {
      p <- qtl2::plot_coef(x            = c2eff,
                           columns      = ncol(c2eff),
                           col          = palette[seq_len(ncol(c2eff))],
                           map          = pmap,
                           scan1_output = scan1_out,
                           bgcolor      = "gray95")
      if (psave) {
        if (is.null(pname)) {
          pname <- paste0(tissue, "_", pheno, "_effects.png")
        }
        if(!is.null(paltte)) {
          png(paste0(outdir, "/", tissue, "_", pheno, "_effects.png"),
              width = 1024, height = 1024)
          print(p)
          dev.off()
        }

        if(!is.null(founders) & pop %notin% c("do", "cc")) {
          legend_plot <- paste0(outdir,"/founder_effects_legend.png")
          if (!file.exists(legend_plot) & is.null(palette)) {
            if (length(founders) <= 3) {
              cols <- c("slateblue","violetred", "green3")
            }
            if (length(founders) > 3 & length(founders) <= 8) {
              cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                        "#66A61E", "#E6AB02", "#A6761D", "#666666")
            }
            png(legend_plot, bg = "transparent")
            plot.new()
            legend("center", founders, lty = 1, col = cols[seq_len(length(founders))],
                   bty = "n", ncol = 4, lwd = 1.3)
            dev.off()
          }
          if (!file.exists(legend_plot) & !is.null(palette)) {
            png(legend_plot, bg = "transparent")
            plot.new()
            legend("center", founders, lty = 1, col =
                     palette[seq_len(length(founders))], bty = "n", ncol = 4,
                   lwd = 1.3)
            dev.off()
          }
        }
      }
    }

  } else {
    p <- plot(scan1_out, pmap)
    if(psave) {
      if(is.null(pname)) {
        pname <- paste0(tissue, "_", pheno, "_LOD.png")
      }
      png(paste0(outdir, "/", pname), width = 1024, height = 730)
      print(p)
      dev.off()
    }
  }
  return(p)
}
