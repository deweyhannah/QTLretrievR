setwd("/projects/munger-lab/projects/eqtl_package/QTLretrievR/")

outdir <- "../RPE_test"
geno_out <- "RPE_gbrs_interpolated_genoprobs.rds"
gbrs_fileLoc <- "/projects/munger-lab/projects/DO_RPE/results/Bowtie_gbrs.out/trim/GBRS_GRCm39/gbrs_tsv_files"
tissues <- c("retina","rpe")
genoprobably(outfile = paste0(outdir,"/",geno_out), gbrsFileLoc = gbrs_fileLoc, tissues = tissues)
