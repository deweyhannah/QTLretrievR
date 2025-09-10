test_that("mapping works correctly", {
  skip_on_cran()
  map_peaks <- mapQTL(genoprobs = list("mESC" = demo_probs),
                      samp_meta = demo_meta,
                      expr_mats = list("mESC" = demo_counts[1:500,]),
                      covar_factors = "sex",
                      gridFile = gridfile69k,
                      annots = demo_annot,
                      save = "ro")

  expect(is.list(map_peaks), "mapping output is not a list")

  expect_in(names(map_peaks), c("maps_list", "peaks_list"))

  expect(is.list(map_peaks$maps_list), "maps_list is not a list")

  expect_in(names(map_peaks$maps_list), c("covar_list", "expr_list",
                                          "exprZ_list", "qtlprobs",
                                          "kinship_loco", "gmap", "pmap",
                                          "map_dat2", "tissue_samp"))

  expect(is.list(map_peaks$peaks_list), "peaks_list is not a list")

  expect(is.data.frame(map_peaks$peaks_list$mESC), "peaks_list does not contain
         dataframes")

  expect_in(colnames(map_peaks$peaks_list$mESC), c("phenotype", "peak_chr",
                                                   "peak_bp", "lod", "ci_lo",
                                                   "ci_hi", "symbol", "chr",
                                                   "start", "end", "midpoint",
                                                   "local"))
})

test_that("error when covariates not in metadata are passed", {
  expect_error(
    mapQTL(genoprobs = list("mESC" = demo_probs),
           samp_meta = demo_meta,
           expr_mats = list("mESC" = demo_counts[1:500,]),
           covar_factors = c("sex","ngen"),
           gridFile = gridfile69k,
           annots = demo_annot,
           save = "ro"),
    "Chosen factors are not in sample metadata. Please check factors and
         sample metadata for missing or misspelled elements"
  )
})

test_that("warning that not all annotations present with wrong annotations", {
  expect_message(
    mapQTL(genoprobs = list("mESC" = demo_probs),
           samp_meta = demo_meta,
           expr_mats = list("mESC" = demo_counts[1:500,]),
           covar_factors = "sex",
           gridFile = gridfile69k,
           annots = annot_105,
           save = "ro"),
    "Not all phenotypes in mESC are present in
                       annotations file. Not all phenotypes will be
                       annotated."
  )
})

test_that("error that not all annotations present with transposed counts", {
  expect_error(
    mapQTL(genoprobs = list("mESC" = demo_probs),
           samp_meta = demo_meta,
           expr_mats = list("mESC" = t(demo_counts[1:500,])),
           covar_factors = "sex",
           gridFile = gridfile69k,
           annots = demo_annot,
           save = "ro"),
    "Phenotypes in mESC not present in annotations file.
                    Please check phenotype names in counts."
  )
})

test_that("error when save with no output directory", {
  expect_error(
    mapQTL(genoprobs = list("mESC" = demo_probs),
           samp_meta = demo_meta,
           expr_mats = list("mESC" = demo_counts[1:500,]),
           covar_factors = "sex",
           gridFile = gridfile69k,
           annots = demo_annot,
           save = "so"),
    "Requested Save. No output directory provided, and no default."
  )
})

test_that("error that not all annotations present with non-rankZ counts", {
  expect_error(
    mapQTL(genoprobs = list("mESC" = demo_probs),
           samp_meta = demo_meta,
           expr_mats = list("mESC" = demo_counts[1:500,]),
           covar_factors = "sex",
           gridFile = gridfile69k,
           annots = demo_annot,
           save = "ro",
           rz = T),
    "Phenotypes in mESC not present in annotations file.
                    Please check phenotype names in counts."
  )
})

test_that("message when rankZ provided and transposed counts", {
  expect_message(
    mapQTL(genoprobs = list("mESC" = demo_probs),
           samp_meta = demo_meta,
           expr_mats = list("mESC" = t(demo_counts[1:500,])),
           covar_factors = "sex",
           gridFile = gridfile69k,
           annots = demo_annot,
           save = "ro",
           rz = T),
    "rankZ transformed counts provided"
  )
})

test_that("annotation columns missing from peak dfs when annotation omitted", {
  map_peaks <- mapQTL(genoprobs = list("mESC" = demo_probs),
                      samp_meta = demo_meta,
                      expr_mats = list("mESC" = demo_counts[1:500,]),
                      covar_factors = "sex",
                      gridFile = gridfile69k,
                      save = "ro")

  expect_false(all(c("phenotype", "peak_chr", "peak_bp","lod", "ci_lo", "ci_hi",
                 "symbol", "chr", "start", "end", "midpoint", "local") %in%
                 colnames(map_peaks$peaks_list$mESC)))
})


test_that("error that not all annotations present with non-rankZ counts", {
  expect_error(
    mapQTL(genoprobs = list("mESC" = demo_probs),
           samp_meta = demo_meta,
           expr_mats = demo_counts[1:500,],
           covar_factors = "sex",
           gridFile = gridfile69k,
           annots = demo_annot,
           save = "ro",
           rz = T),
    "More expression matrices than probabilites found. Are you using the
         correct probabilities?"
  )
})
