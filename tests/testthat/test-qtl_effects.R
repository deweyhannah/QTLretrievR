test_that("effects calling works with correct inputs", {
  effects_out <- qtl_effects(peaks = map_peaks_obj$peaks_list,
                             mapping = map_peaks_obj$maps_list,
                             save = "ro")

  expect(is.list(effects_out), "effects output not a list")
  expect_named(effects_out, c("effects_blup", "peaks"))

  expect_equal(colnames(effects_out$effects_blup$mESC), LETTERS[1:8])
  expect_equal(colnames(effects_out$peaks$mESC), c("phenotype", "peak_chr",
                                              "peak_bp", "lod", "ci_lo",
                                              "ci_hi", "symbol", "chr",
                                              "start", "end", "midpoint",
                                              "local", "before", "after",
                                              "marker"))
})
