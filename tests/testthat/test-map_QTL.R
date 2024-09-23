
message( timestamp())
ptm <- proc.time()
outdir <- ("/projects/munger-lab/projects/SA_eqtl_package/RPE_test")

genoprobs <- "retina_genoprobs.rds"#"RPE_gbrs_interpolated_genoprobs_onejob.rds"
peaks_out <- "RPE_mm39_peaks.rds"
mapping_out <- "RPE_mm39_mapping.rds"
metadata <-  "/projects/munger-lab/projects/eqtl_package/RPE_test/data/sample_table.txt"
expr_mats <- c("/projects/munger-lab/projects/SA_eqtl_package/RPE_test/retina_normalized_counts_subset.txt")
               #,"/projects/munger-lab/projects/SA_eqtl_package/RPE_test/RPE_normalized_counts_subset.txt")
# expr_mats <- c("/projects/munger-lab/projects/eqtl_package/RPE_test/data/retina_normalized_counts.txt", "/projects/munger-lab/projects/eqtl_package/RPE_test/data/RPE_normalized_counts.txt")
covar_factors <- c("sex", "diet")
biomart <- "/projects/munger-lab/projects/DO_RPE/ensembl_v105_mouse_genes.txt"


map_peaks <- mapQTL(outdir, peaks_out, mapping_out, genoprobs, metadata, expr_mats, covar_factors, biomart = biomart, n.cores = 8)

message( timestamp())
message(proc.time() - ptm)
message(paste0(names(map_peaks), sep = " "))
message(paste0(names(map_peaks$maps_list), sep = " "))
message(paste0(names(map_peaks$peaks_list), sep = " "))
