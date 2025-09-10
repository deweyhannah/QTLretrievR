try({
  map_peaks_obj <- mapQTL(genoprobs = list("mESC" = demo_probs),
                          samp_meta = demo_meta,
                          expr_mats = list("mESC" = demo_counts[1:500,]),
                          covar_factors = "sex",
                          gridFile = gridfile69k,
                          annots = demo_annot,
                          save = "ro")
}, silent = TRUE)

try({
  med_obj <- modiFinder(peaks = map_peaks_obj$peaks_list,
                        mapping = map_peaks_obj$maps_list,
                        annots = demo_annot,
                        save = "ro")
}, silent = TRUE)
