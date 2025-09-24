---
title: "QTL Identification and Analysis with QTLretrievR"
output: 
  rmarkdown::html_vignette:
    toc: true
    toc_depth: 2
vignette: >
  %\VignetteIndexEntry{QTL Identification and Analysis with QTLretrievR}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

<div align='center'>**Abstract**</div>
Advances in molecular phenotyping technologies, including quantitative transcriptomics and proteomics (“-omics”), have transformed biomedical research over the past two decades by providing a granular understanding of the gene regulatory networks that drive disease. Many studies have demonstrated the additional power of combining  -omics profiling with genetic mapping, an integrative approach referred to as “systems genetics”, to link population variability in complex disease phenotypes back to genetic variants and their direct effects on gene regulation (e.g. gene expression). However, these quantitative trait locus (QTL) analyses are often computationally intensive and intimidating to those without a strong computational background — especially when samples originate from multiple cell/tissue types or involve thousands of molecular phenotypes. We have developed a new R package, QTLretrievR, that aims to simplify and streamline the identification and downstream analyses of molecular QTL from genetically diverse populations such as Diversity Outbred or Collaborative Cross mice. QTLretrievR leverages existing R packages used in complex trait analysis, including the popular `qtl2` package for QTL peak calling and haplotype effect identification and the `intermediate` package for peak mediation, and employs dynamic nested parallelism to enhance the efficiency of each analysis step. Future developments will enhance the visualization functionality and data analysis options, including adding network inference methods and expanding the software to work with other mapping populations.




```r
library(QTLretrievR)
```

# Quick start

This is the easiest way to run QTLretrievR. There are several input files required for `runQTL` as described in detail below. This function will return the list of all the objects used in QTL mapping, annotated QTL peaks, QTL effects, and mediation results. Note in the example below that the objects will be only be returned, not saved. If you wish to have the objects saved for later use you can change `save` to "so" (save only), or "sr" (save and return).


```r
## Passing a Genoprobs object - list of probabilities for each tissue
## Using internal data for expr, gridFile, metadata, annots
## Unevaluated Code Chunk

qtl_res <- runQTL(outdir = "../../vignette/", 
                  gbrs_fileLoc = list("mESC" = demo_probs),
                  geno_out = "mESC_mm10_probs.rds", 
                  peaks_out = "mESC_mm10_peaks.rds", 
                  map_out = "mESC_mm10_map.rds", 
                  med_out = "mESC_mm10_mediation.rds",
                  effects_out = "mESC_mm10_effects.rds",
                  expr_mats = list("mESC" = demo_counts),
                  tissues = "mESC", 
                  gridFile = gridfile69k, 
                  covar_factors = c("sex"),
                  metadata = demo_meta, 
                  annots = demo_annot, 
                  save_t = "ro")
```

This can be alternatively submitted as follows:


```r
## Unevaluated Code Chunk

qtl_res <- runQTL(outdir = <path/to/output/dir>, 
                  gbrs_fileLoc = <path/to/gbrs/files>,
                  geno_out = "Name_for_genotype_prob_file.rds", 
                  peaks_out = "Name_for_QTL_peaks.rds", 
                  map_out = "Name_for_mapping_objects.rds", 
                  med_out = "Name_for_mediation_Results.rds",
                  effects_out = "Name_for_QTL_effects.rds",
                  expr_mats = c(<path/to/expr_mat1>, <path/to/expr_mat2>),
                  tissues = "Tissue_name",
                  gridFile = <path/to/gridfile>, 
                  covar_factors = c("Covariates"),
                  metadata = <path/to/sample/metadata>, 
                  annots = <path/to/annotations>, 
                  save = "ro"
                  )

```

This returns the mapping object, original peaks object, mediation object, and effects object. The peaks, effects and mediation objects can be used for additional plotting and analysis.

# Inputs Data

There are several types of data that `QTLretrievR` takes as input, below we describe variations on the required and optional data inputs.

## Genotype Probabilities

Genotype Probabilities (genoprobs) can be passed to `QTLretreivR` in a few different ways. 

  - Passed as a string to the location of the .tsv probability files from the GBRS software. This can be the folder where all the GBRS results are kept, the files should end with `.gbrs.interpolated.genoprobs.tsv`. This option would be used with `genoprobably` and `runQTL`. This is probably the easiest option if you are starting from scratch and only have GBRS genotypes.
  
  - Passed as a two dimensional array of the genotype probabilities for all the samples. See the `genotype_probs.Rds` file [here](https://figshare.com/articles/dataset/Supplemental_data_repository_for_Skelly_et_al_2020_Cell_Stem_Cell_Mapping_the_effects_of_genetic_variation_on_chromatin_state_and_gene_expression_reveals_loci_that_control_ground_state_pluripotency_/12233570) for an example of how this should look. This option is for passing to `runQTL`.
  
  - Passed as a `qtl2` formatted genoprobs object. This option is for passing to `runQTL` and `mapQTL`. An example of this is shown below.
  

```r
## Structure and example of qtl2 fomatted genoprobs
names(demo_probs)
#>  [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "X"
str(demo_probs[[1]])
#>  num [1:10, 1:8, 1:4711] 1.48e-08 1.23e-08 1.04e-08 3.31e-08 2.94e-07 ...
#>  - attr(*, "dimnames")=List of 3
#>   ..$ : chr [1:10] "PB357.02_repA" "PB357.12_repA" "PB357.16_repA" "PB357.18_repA" ...
#>   ..$ : chr [1:8] "A" "B" "C" "D" ...
#>   ..$ : chr [1:4711] "1_3000000" "1_3041392" "1_3346528" "1_3651663" ...
demo_probs[[1]][1:2,,1:2]
#> , , 1_3000000
#> 
#>                          A   B            C            D            E            F            G            H
#> PB357.02_repA 1.477899e-08 0.5 1.217364e-08 1.370868e-08 1.199720e-08 5.000000e-01 4.331041e-09 1.313560e-08
#> PB357.12_repA 1.232415e-08 0.5 5.000000e-01 1.604246e-08 3.499095e-09 1.478095e-08 4.742920e-09 1.118048e-08
#> 
#> , , 1_3041392
#> 
#>                          A   B            C            D            E            F            G            H
#> PB357.02_repA 1.477899e-08 0.5 1.217364e-08 1.370868e-08 1.199720e-08 5.000000e-01 4.331041e-09 1.313560e-08
#> PB357.12_repA 1.232415e-08 0.5 5.000000e-01 1.604246e-08 3.499095e-09 1.478095e-08 4.742920e-09 1.118048e-08
```
  
If you have genotyping from one of the Mouse Universal Genotyping Arrays (MUGA), then you can format your probabilities using our `mugaprobs` function, a detailed example is below in the [Genoprobs](#genoprobs) section of this document.
  
## Quantified Phenotypes

Although the variable is named `expr_mats` any quantified molecular data can be passed to this variable. They can be passed as a vector of paths to the data (in the same order as tissues are being passed), or as a named list of matrices (named with the tissue type associated with the counts). These should be filtered, batch corrected and normalized data with the phenotype in rows and the samples in columns, there is no need to rankZ transform your counts, that step occurs in `mapQTL`. See the attached `demo_counts` for an example of formatting.


```r
## Information and partial example of quantified data
class(demo_counts)
#> [1] "matrix" "array"
head(demo_counts)
#>                    PB357.02_repA PB357.12_repA PB357.16_repA PB357.18_repA PB357.19_repA PB357.20_repA PB357.21_repA
#> ENSMUSG00000030402      90.27620      47.49881      65.92215      60.68472       50.6948      43.98085      68.98128
#> ENSMUSG00000025473      55.33058      61.24847      59.68627      26.61611       24.3725      45.02802      55.41119
#> ENSMUSG00000029246    2333.74400    5135.62454    2341.10901    2701.05069     3481.3853    2611.54250    1373.65252
#> ENSMUSG00000051413     218.41017     162.49593     209.34736     223.57530      230.0764     158.12164     211.46720
#> ENSMUSG00000026356    3244.11907   11302.21694    2892.56671    4801.55973     5616.3992    4653.59323    2312.56906
#> ENSMUSG00000042505     317.42278     554.98610     316.24815     300.22969      348.0393     347.65818     360.73816
#>                    PB357.22_repA PB357.28_repA PB357.29_repA
#> ENSMUSG00000030402      57.27735      42.11575      52.36949
#> ENSMUSG00000025473      37.55891      14.03858      44.38093
#> ENSMUSG00000029246    3451.91477    3130.74325    2534.28661
#> ENSMUSG00000051413     217.84171     161.91166     260.95984
#> ENSMUSG00000026356    5754.02576    6941.61141    4452.29448
#> ENSMUSG00000042505     294.83748     440.81151     286.70078
```

## Marker Grid

This should be a grid of markers and the associated genomic locations in bp and cM. The grid file can be passed as an object, or as a path to where the file is located. There are two grid files attached to `QTLretrievR`, a grid with ~75K pseudo-markers with positions updated for mm38 and a grid with ~69K pseudo-markers with positions associated with mm10, and the GigaMUGA and MegaMUGA markers and positions. Or you can use a  new marker grid to be used with your probabilities. If making your own grid file, see `gridfile` to make sure you have the correct format.


```r
## example format for 75K marker grid
head(gridfile)
#>           chr     pos      cM
#> 1_3000000   1 3000000 0.00001
#> 1_3039563   1 3039563 0.02001
#> 1_3079126   1 3079126 0.04001
#> 1_3118689   1 3118689 0.06001
#> 1_3158252   1 3158252 0.08001
#> 1_3197814   1 3197814 0.10001
```

## Metadata

At a minimum, the metadata file should contain a column for sample IDs (`ID`) and columns for any covariates to be included in the peak and effect calling. Any other pieces of information that may be relevant to your analysis that you want to make sure is all in the same place should also be included. A minimal example, `demo_meta`, is attached to `QTLretrievR`.


```r
## minimal metadata example
head(demo_meta[which(demo_meta$ID %in% colnames(demo_counts)),])
#>              ID sex
#> 1 PB357.02_repA   F
#> 2 PB357.12_repA   F
#> 3 PB357.16_repA   F
#> 4 PB357.18_repA   F
#> 5 PB357.19_repA   F
#> 6 PB357.20_repA   F
```

## Annotations

Currently the annotations needs to contain at a minimum `id`, `symbol`, `chr`, `start`, and `end` columns. The `id` column should contain the phenotype information that appears in the row names of the quantified data. In the attached example (`demo_annot`) the id column contains ensembl gene IDs, these could be traded out for ensembl protein IDs or any other phenotype identifier. The `demo_annot` file is the Ensemblv84 biomart annotations, if you are using newer expression data use `annot_105` or upload your own annotations


```r
## Example annotation file - Note this file is all the genes for the Ensemblv84 build
head(demo_annot)
#>                   id symbol chr     start       end
#> 1 ENSMUSG00000000001  Gnai3   3 108107280 108146146
#> 2 ENSMUSG00000000028  Cdc45  16  18780447  18811987
#> 3 ENSMUSG00000000031    H19   7 142575529 142578143
#> 4 ENSMUSG00000000037  Scml2   X 161117193 161258213
#> 5 ENSMUSG00000000056   Narf  11 121237253 121255856
#> 6 ENSMUSG00000000058   Cav2   6  17281185  17289115
```

# Running as Separate Steps

If, instead of using the wrapper to run all steps at once, you wish to run each step individually (calculating/formatting the genotype probabilities, mapping/peak calling, mediation, effects) this can be done with relative ease.

## Genoprobs {#genoprobs}

In order to convert the genoprobs from the GBRS output to the `qtl2` format we run a function called `genoprobably`. It pulls in the relevant files from the provided location and splits samples into the tissues they are derived from (the names that you pass to `tissue` should be included in the file name in some form), then puts these into a 3 dimensional format based on the markers. Once the genotype probabilities are in the 3D format `probs_3d_to_qtl2` is used to convert them to the qtl2 format. The below example shows how to run genoprobably with the 75k gridfile.


```r
## Unevaluated Code Chunk
probs <- genoprobably(outfile     = "../../vignette/gbrs_interpolated_probs.rds", 
                      gbrsFileLoc = "<path/to/gbrs/results>",
                      tissues     = c("tissue1","tissue2"), 
                      gridFile    = gridfile, 
                      save        = "ro")
```

Alternatively you can use genotyping from either GigaMUGA or MegaMUGA as your genome probabilities. In order to convert these reports to the `qtl2` format we run a function called `mugaprobs`. It loads the appropriate grid files, builds the control file, and processes the genotyping file(s) into the correct allele based format for processing. The below example shows how to run `mugaprobs` from a sample GigaMUGA genotyping run.


```r
## Unevaluated Code Chunk

probs <- mugaprobs(type       = "GM",                                   # Did your genotyping come from GigaMUGA or MegaMUGA?
                   covarLoc   = "/path/to/covariate/file/", 
                   covar_file = "name_of_covariate_file.tsv", 
                   i.files    = c("final_report_1", "final_report_2"),  # If you already have chromosome specific genotype files you can pass that directory as a string here
                   genoPrefix = "gm4qtl2", 
                   probsOut   = "muga_interpolated_genoprobs.rds",      # This is the file that your probabilites will be saved to
                   tissues    = c("tissue1", "tissue2"))                # All probabilities will be attributed to each tissue, this is just to put it in the format needed for `mapQTL`

```

## Mapping and Peak Calling

The mapping and peak calling steps are combined in our pipeline. We use the genoprobs to calculate kinship, and create the genetic and pysical maps. Then using the quantified data and sample metadata we rankZ normalize the counts, calculate covariates and finally do the peak calling using `qtl2::scan1`.


```r
## Unevaluated Code Chunk
## Passing objects:
map_peaks <- mapQTL(outdir        = "../../vignette/", 
                    peaks_out     = "mm39_peaks.rds",
                    map_out       = "mm39_map.rds", 
                    genoprobs     = probs, 
                    samp_meta     = demo_meta,                   # Make sure that your ID column is named "ID"
                    expr_mats     = list("mESC" = demo_counts),  # Pass your counts as a named list
                    covar_factors = c("sex"),
                    gridFile      = gridfile, 
                    annots        = demo_annot,                  # See the Annotations section above for notes on this
                    save          = "ro")                        # return only - other options are "so" (save only) and "sr" (save and return)

## Passing file locations
map_peaks <- mapQTL(outdir        = "../../vignette/", 
                    peaks_out     = "mm39_peaks.rds",
                    map_out       = "mm39_map.rds", 
                    genoprobs     = "path/to/saved/probs.rds",
                    samp_meta     = "path/to/sample/metadata",
                    expr_mats     = c("path/to/counts1","path/to/coounts/2"), 
                    covar_factors = c("sex"),
                    gridFile      = "path/to/custom/grid",
                    annots        = "path/to/annotations",
                    save          = "ro")
```

## Mediation

For mediation we use the peaks, mapping information (which includes the rankZ normalized counts), and annotations to identify targets which in turn allows us to identify phenotypes that act as mediators for peaks above a suggestive LOD. This is done using `intermediate::mediation.scan`.


```r
## Unevaluated Code Chunk
med_res <- run_mediate(peaks   = map_peaks$peaks_list, 
                       mapping = map_peaks$maps_list,
                       outdir  = "../../vignette", 
                       annots  = demo_annot,
                       suggLOD = 7, 
                       med_out = "mm39_mediation.rds", 
                       save    = "ro")
```

## Effects

To identify the founder effects at peaks of interest, we use the peaks, and mapping information along with a suggestive LOD threshold to ensure that the reported effects are for relevant peaks only. This is done using `qtl2::scan1blup`.


```r
## Unevaluated Code Chunk
effects <- qtl_effects(mapping = map_peaks$maps_list, 
                       peaks   = map_peaks$peaks_list,
                       suggLOD = 7, 
                       outdir  = "../../vignette", 
                       outfile = "mm39_effects.rds", 
                       save    = "ro")
```


# Submitting to a Job Scheduler

Currently, `QTLretrievR` is set up to identify available cores in a linux, mac, or windows environment. It is also able to identify the number of cores if you are submitting a job through a SLURM scheduler. Just make sure to include both of the following lines in your `SBATCH`:

```{}
#SBATCH --ntasks=N
#SBATCH --cpus-per-task=X
```

When the number of cores is calculated for parallelization purposes, the total number of cores that are identified are `N*X`. If you don't include both, and are submitting with a SLURM scheduler, it may not identify the correct number of cores.

If you are using any other scheduler, make sure that you specify the total number of cpus for the whole job. If applicable the number of processes should be 1.
