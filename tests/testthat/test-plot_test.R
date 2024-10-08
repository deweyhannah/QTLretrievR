setwd("/projects/munger-lab/projects/eqtl_package/QTLretrievR/")

## Output directory
outputdir <- "../RPE_test"

## Read in suggestive peaks
effects <- readRDS("../RPE_test/RPE_mm39_effects_onejob.rds")

## Read in map data
mapping <- readRDS("../RPE_test/RPE_mm39_mapping_onejob.rds")
map_dat2 <- mapping$map_dat2

## Clean up
rm(mapping)

qtl_map <- plot_eqtlmap(map_dat = map_dat2, peaks = effects$peaks, outdir = outputdir, outbase = "RPE_mm39_eqtl_map", unit = "mbp")

hotspots <- transbands(map_dat = map_dat2, peaks = effects$peaks, outdir = outputdir)

message(paste0(names(qtl_map), sep = " "))
message(paste0(names(hotspots), sep = " "))

message(paste0(unique(hotspots$bands.rna$rpe$chr), sep = " "))
message(paste0(dim(hotspots$bands.rna$rpe), sep = " "))
