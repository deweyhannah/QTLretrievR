message( timestamp())
ptm <- proc.time()
## Output directory
outputdir <- ("/projects/munger-lab/projects/SA_eqtl_package/RPE_test")

## Output files
geno <- "RPE_gbrs_interpolated_genoprobs_onejob.rds"
peaks <- "RPE_mm39_peaks_onejob.rds"
mapping <- "RPE_mm39_mapping_onejob.rds"
mediation <- "RPE_mm39_mediation_onejob.rds"
effects <- "RPE_mm39_effects_onejob.rds"

## Input Files
gbrs_files <- "/projects/munger-lab/projects/DO_RPE/results/Bowtie_gbrs.out/trim/GBRS_GRCm39/gbrs_tsv_files"
sample_data <- "/projects/munger-lab/projects/eqtl_package/RPE_test/data/sample_table.txt"
exprs <-  c("/projects/munger-lab/projects/SA_eqtl_package/RPE_test/retina_normalized_counts_subset.txt","/projects/munger-lab/projects/SA_eqtl_package/RPE_test/RPE_normalized_counts_subset.txt")
# c("/projects/munger-lab/projects/eqtl_package/RPE_test/data/retina_normalized_counts.txt", "/projects/munger-lab/projects/eqtl_package/RPE_test/data/RPE_normalized_counts.txt")
biomart_file <- "/projects/munger-lab/projects/DO_RPE/ensembl_v105_mouse_genes.txt"

## Other info
covars <- c("sex", "diet")
tissues_list <- c("retina","rpe")

ret_excl <- read.delim("/projects/munger-lab/projects/eqtl_package/RPE_test/data/retina_excl_samples.txt", skip = 1, header = F, sep = " ")
rpe_excl <- read.delim("/projects/munger-lab/projects/eqtl_package/RPE_test/data/rpe_excl_samples.txt", skip = 1, header = F, sep = " ")
exclude <- c(ret_excl$V2, rpe_excl$V2)


qtl_all <- runQTL(geno_out = geno,
                  gbrs_fileLoc = gbrs_files,
                  tissues = tissues_list,
                  outdir = outputdir,
                  peaks_out = peaks,
                  map_out = mapping,
                  med_out = mediation,
                  effects_out = effects,
                  metadata = sample_data,
                  expr_mats = exprs,
                  covar_factors = covars,
                  n.cores = 8,
                  suggLOD = 7,
                  biomart = biomart_file,
                  localRange = 10e6
                  )
message( timestamp())
message(proc.time() - ptm)
message(paste0(names(qtl_all), sep = " "))
message(paste0(names(qtl_all$peaks_list), sep = " "))
message(paste0(names(qtl_all$maps_list), sep = " "))
message(paste0(names(qtl_all$res_list), sep = " "))
message(paste0(names(qtl_all$effects_res), sep = " "))


