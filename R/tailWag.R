#' Create a plot with founder effects for a given phenotype
#'
#' @param effects Effects object (output from `qtl_effects` function).
#' @param tissue String indicating which tissue to plot the founder effects for.
#' @param feat Individual phenotype to plot effects of.
#' @param chromosome Chromosome that the peak is present on.
#' @param color Plot color. Default is "#0073C2FF".
#' @param psave Logical. Save the plot as `.png`. Default `TRUE`.
#' @param pname File name to save plot (needs to end in `.png`).
#'  Default is `wag_effects_<tissue>_<feat>.png`
#' @param outdir Directory to save plots. Default is `NULL`.
#' @param pop One of `c("do", "cc", "other")` to indicate founder population.
#'  Default is "`do`".
#' @param founders If `pop == "other"`, list of founders in haplotype order.
#' @param symb Logical. The phenotype to be plotted is passed as a symbol
#' instead of phenotype id. Default `TRUE`.
#'
#' @return ggplot object of the effects plot
#'
#' @export
#'
tailWag <- function(effects, tissue, feat, chromosome, color = "#0073C2FF",
                    psave = TRUE, pname = NULL, outdir = NULL, pop = "do",
                    founders = NULL, symb = TRUE) {
  if (psave & is.null(outdir)) {
    stop("Plot to be saved, but no directory provided")
  }
  if (psave & !is.null(pname)) {
    pname <- paste0("wag_efects_", tissue, "_", feat, ".png")
    message(paste0("Plot to be saved. Saving as ", pname, " in ", outdir))
  }

  if (pop %in% c("do", "cc")) {
    founders <- c("AJ","B6","129","NOD","NZO","CAST","PWK","WSB")
  }

  peak_effects <- cbind(effects$peaks[[tissue]], effects$effects_blup[[tissue]])

  # message(colnames(peak_effects), sep = "\t")

  # if (symb) {
  #   feat_effects <- peak_effects |>
  #     dplyr::filter(symbol %in% feat,
  #                   peak_chr == chromosome) |>
  #     dplyr::rename(feature = symbol) |>
  #     dplyr::select(feature, LETTERS[seq_len(length(founders))])
  # } else {
    feat_effects <- peak_effects |>
      dplyr::filter(phenotype %in% feat,
                    peak_chr == chromosome) |>
      dplyr::rename(feature = phenotype) |>
      dplyr::select(feature, LETTERS[seq_len(length(founders))])
  # }

  colnames(feat_effects) <- c("feature", founders)

  feat_effects_df <- feat_effects  |>
    tidyr::pivot_longer(cols      = founders,
                        names_to  = "effect",
                        values_to = "value") |>
    dplyr::mutate(feature = factor(feature, levels = feat),
                  effect = factor(effect, levels = rev(founders))) |>
    as.data.frame()
  ymin <- floor(min(feat_effects_df$value))
  ymax <- ceiling(max(feat_effects_df$value))

  n_feats <- length(unique(feat_effects_df$feature))
  if (length(color) == 1 && n_feats > 1) {
    color = c(color, shift_hue(color))
    names(color) <- feat
  }

  p <- feat_effects_df |>
    ggplot2::ggplot() +
    ggplot2::aes(x = effect,
                 y = value,
                 color = feature,
                 group = feature) +
    ggplot2::geom_point(size = 4, show.legend = FALSE) +
    ggplot2::geom_line(size = 1.2, show.legend = TRUE) +
    ggplot2::theme_minimal(base_size = 18) +
    ggplot2::ylab("Haplotype effects") +
    ggplot2::xlab("") +
    ggplot2::ylim(ymin, ymax) +
    ggplot2::scale_color_manual(values = color, name = "Phenotype") + #setNames(color, feat)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::theme(axis.line.x = element_blank(),
                   axis.title = element_text(size = 18)) +
    ggplot2::coord_flip( clip ="off") +
    ggplot2::theme(legend.position = "top")

  if (psave) {
    ggplot2::ggsave(pname, plot = p, device = "png", path = outdir,
                    bg = "white")
  }
  return(p)
}

shift_hue <- function(hex_color, shift = 180) {
  rgb_vals <- grDevices::col2rgb(hex_color) / 255

  # RGB to HSL
  cmax <- max(rgb_vals); cmin <- min(rgb_vals); delta <- cmax - cmin
  l <- (cmax + cmin) / 2
  s <- if (delta == 0) 0 else delta / (1 - abs(2 * l - 1))
  h <- if (delta == 0) 0 else switch(
    which.max(rgb_vals),
    `1` = 60 * (((rgb_vals[2] - rgb_vals[3]) / delta) %% 6),
    `2` = 60 * (((rgb_vals[3] - rgb_vals[1]) / delta) + 2),
    `3` = 60 * (((rgb_vals[1] - rgb_vals[2]) / delta) + 4)
  )

  # Shift hue and convert back to RGB
  h <- (h + shift) %% 360
  c_val <- (1 - abs(2 * l - 1)) * s
  x <- c_val * (1 - abs((h / 60) %% 2 - 1))
  m <- l - c_val / 2
  rgb_new <- switch(
    floor(h / 60) + 1,
    c(c_val, x, 0), c(x, c_val, 0), c(0, c_val, x),
    c(0, x, c_val), c(x, 0, c_val), c(c_val, 0, x)
  )
  grDevices::rgb(rgb_new[1] + m, rgb_new[2] + m, rgb_new[3] + m)
}
