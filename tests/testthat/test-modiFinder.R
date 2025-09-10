test_that("mediation works with correct inputs", {
  demo_med <- modiFinder(peaks = map_peaks_obj$peaks_list,
                         mapping = map_peaks_obj$maps_list,
                         annots = demo_annot,
                         save = "ro")

  expect(is.list(demo_med), "mediation output is not a list")
  expect_equal(names(demo_med), "mESC")

  expect_setequal(colnames(demo_med$mESC), c("target_id", "qtl_lod", "qtl_chr",
                                             "mediator", "mediator_id",
                                             "mediator_chr",
                                             "mediator_midpoint", "LOD",
                                             "target"))
})

test_that("hotspot only mediation works with correct inputs", {
  demo_med_hs <- modiFinder(peaks = map_peaks_obj$peaks_list,
                            mapping = map_peaks_obj$maps_list,
                            annots = demo_annot,
                            save = "ro",
                            hsOnly = T)

  expect(is.list(demo_med_hs), "mediation output is not a list")
  expect_equal(names(demo_med_hs), "mESC")

  expect_setequal(colnames(demo_med_hs$mESC), c("target_id", "qtl_lod",
                                                "qtl_chr","mediator",
                                                "mediator_id", "mediator_chr",
                                                "mediator_midpoint", "LOD",
                                                "target"))

  expect_lt(nrow(demo_med_hs$mESC), nrow(med_obj$mESC))
})
