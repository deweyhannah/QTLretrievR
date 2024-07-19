setwd("/projects/munger-lab/projects/eqtl_package/QTLretrievR/")

outdir <- "../RPE_test"
genoprobs <- "RPE_gbrs_interpolated_genoprobs.rds"
peaks_out <- "RPE_mm39_peaks.rds"
mapping_out <- "RPE_mm39_mapping.rds"
metadata <- "../RPE_test/data/sample_table.txt"
expr_mats <- c("../RPE_test/data/retina_normalized_counts.txt", "../RPE_test/data/RPE_normalized_counts.txt")
covar_factors <- c("sex", "diet")
biomart <- "/projects/munger-lab/projects/DO_RPE/ensembl_v105_mouse_genes.txt"

map_peaks <- mapeQTL(outdir, peaks_out, map_out, genoprobs, metadata, expr_mats, covar_factors, biomart = biomart, samp_excl = c())
list2env(map_peaks)
ls()
