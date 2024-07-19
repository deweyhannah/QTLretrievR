## usethis namespace: start
#' @export
## usethis namespace: end
runQTL <- function(geno_out = "gbrs_interpolated_genoprobs.RDS", gbrs_fileLoc, tissues = c(), gridfile = "/projects/compsci/omics_share/mouse/GRCm39/supporting_files/emase_gbrs/rel_2112_v8/ref.genome_grid.GRCm39.tsv",
                   outdir, peaks_out = "mm39_peaks.RDS", map_out = "mm39_mapping.RDS", med_out = "mm39_mediation_res.RDS", effects_out = "mm39_effects.RDS", metadata, expr_mats, covar_factors, n.cores = 4, suggLOD = 7, biomart, localRange = 10e6,
                   samp_excl = c()){

  ## Check oudir
  if(length(outdir) == 0 | !dir.exists(outdir)){
    message("Invalid or no directory provided. Making an output file directory in the current working directory.")
    dir.create("./QTL_mapping")
    outdir <- "./QTL_mapping"
  }

  ## Convert genoprobs
  genoprobs <- genoprobably(outfile = paste0(outdir,"/",geno_out), gbrsFileLoc = gbrs_fileLoc, tissues = tissues, gridfile = gridfile)

  ## Map QTLs from genoprobs
  map_peaks <- mapQTL(outdir, peaks_out, map_out, genoprobs, metadata, expr_mats, covar_factors, n.cores, gridfile, localRange, biomart, samp_excl)
  list2env(map_peaks,.GlobalEnv)
  rm(map_peaks)

  ## Run Mediation and Effects
  res_list <- prep_mediate(peaks_list, maps_list, suggLOD, outdir, biomart, med_out)
  effects_res <- qtl_effects(maps_list, peaks_list, suggLOD, outdir, effects_out, n.cores)

  ## For now we are just going to return the peaks, mapping data, mediation results, and effects results
  return(list(peaks_list, maps_list, res_list, effects_res))
}
