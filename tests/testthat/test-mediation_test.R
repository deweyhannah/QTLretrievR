setwd("/projects/munger-lab/projects/eqtl_package/QTLretrievR/")

outdir <- "../RPE_test"
peaks_out <- "RPE_mm39_peaks.rds"
mapping_out <- "RPE_mm39_mapping.rds"
biomart <- "/projects/munger-lab/projects/DO_RPE/ensembl_v105_mouse_genes.txt"
med_out <- "RPE_mm39_mediation.rds"


res_list <- prep_mediate(peaks_out, mapping_out, 7, outdir, biomart, med_out, n.cores = 4)
