**Abstract**

Quantitative trait locus (QTL) analysis is a frequently used technique
in omics studies to identify molecular phenotypes that have effects on
what is being measured. These analyses are often computationally
intensive and intimidating to those without a strong computational
background, especially when samples are taken for multiple cell types
that won’t be analyzed concurrently. QTLretirevR is designed to simplify
and streamline the process of identifying QTL peaks, mediators of those
peaks, and specific local effects. Built on `qtl2` (for peak calling and
effect identification) and `intermediate` (for peak mediation),
QTLretrievR is designed to use nested parallelism in order to optimize
the amount of time it takes to run each analysis step. In its first
version it is designed to work with raw genotype probabilities from
genome reconstruction by RNA-Seq (GBRS) (or any qtl2 formatted genome
probabilities) for all samples, and normalized quantification data from
Collaborative Cross and Diversity Outbred mice.

``` r
library(QTLretrievR)
```

# Standard workflow

## Quick start

This is the easiest way to run QTLretrievR. There are several input
files required for `runQTL` as described in detail below. This function
will return the list of all the objects used in QTL mapping, annotated
QTL peaks, QTL effects, and mediation results. Note in the example below
that the objects will be only be returned, not saved. If you wish to
have the objects saved for later use you can change `save` to “so” (save
only), or “sr” (save and return).

``` r
## Passing a Genoprobs object - list of probabilities for each tissue
## Using internal data for expr, gridFile, metadata, biomart

qtl_res <- runQTL(outdir = "../vignettes/", 
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
                  biomart = demo_annot, 
                  save_t = "ro")
#> using provided probabilities
#> running mapping
#> mESC
#> gmap complete
#> pmap complete
#> expression loaded
#> Working with 10 samples and 1000 genes in mESC.
#> rankZ normalized
#> covariates calculated
#> calculating peaks
#> Working in Linux and there are 40 cores in total.
#> Registering 40 cores and passing 40 cores per tissue to 1 tissue(s). Does that look right? If not please set total_cores parameter to the number of available cores.
#> adding annotations to peaks
#> running mediation
#> load annotations
#> checking peaks and mapping
#> 'data.frame':    19005 obs. of  13 variables:
#>  $ phenotype     : chr  "ENSMUSG00000035984" "ENSMUSG00000015134" "ENSMUSG00000042211" "ENSMUSG00000039089" ...
#>  $ peak_chr      : Factor w/ 20 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ peak_cM       : num  1.46 1.46 1.46 1.46 1.46 ...
#>  $ lod           : num  7.94 9.29 10.19 8.82 9.22 ...
#>  $ ci_lo         : num  1.45 1.45 1.45 1.45 1.45 ...
#>  $ ci_hi         : num  23.2 72.4 21.4 53.1 14.1 ...
#>  $ interp_bp_peak: num  3020696 3020696 3020696 3020696 3020696 ...
#>  $ symbol        : chr  "Nme5" "Aldh1a3" "Fbxo38" "L3mbtl3" ...
#>  $ chr           : chr  "18" "7" "18" "10" ...
#>  $ start         : int  34562634 66390890 62504069 26274468 160645888 123214780 148039077 75175942 65651726 127063534 ...
#>  $ end           : int  34579115 66427517 62548743 26375185 160722764 123273253 148059551 75178835 65772752 127067920 ...
#>  $ midpoint      : num  3.46e+07 6.64e+07 6.25e+07 2.63e+07 1.61e+08 ...
#>  $ local         : num  0 0 0 0 0 0 0 0 0 0 ...
#> data checked
#> filtered peaks
#> running mediation
#> Working in Linux and there are 40 cores in total.
#> Registering 40 cores and passing 40 cores per tissue to 1 tissue(s).
#> running effects
#> data checked
#> peaks extracted, calculating effects now
#> Working in Linux and there are 40 cores in total.
#> Registering 40 cores and passing 40 cores per tissue to 1 tissue(s).
#> effects calculated. saving to RDS
```

This can be alternatively submitted as follows:

``` r
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
                  biomart = <path/to/annotations>, 
                  save = "ro"
                  )
```

This returns the mapping object, original peaks object, mediation
object, and effects object. The peaks, effects and mediation objects can
be used for additional plotting and analysis.

## Input Data

There are several types of data that `QTLretrievR` takes as input, below
we describe variations on the required and optional data inputs.

### Genotype Probabilities

Genotype Probabilities (genoprobs) can be passed to `QTLretreivR` in a
few different ways.

-   Passed as a string to the location of the .tsv probability files
    from the GBRS software. This can be the folder where all the GBRS
    results are kept, the files should end with
    `.gbrs.interpolated.genoprobs.tsv`. This option would be used with
    `genoprobably` and `runQTL`. This is probably the easiest option if
    you are starting from scratch and only have GBRS genotypes.

-   Passed as a two dimensional array of the genotype probabilities for
    all the samples. See the `genotype_probs.Rds` file
    [here](https://figshare.com/articles/dataset/Supplemental_data_repository_for_Skelly_et_al_2020_Cell_Stem_Cell_Mapping_the_effects_of_genetic_variation_on_chromatin_state_and_gene_expression_reveals_loci_that_control_ground_state_pluripotency_/12233570)
    for an example of how this should look. This option is for passing
    to `runQTL`.

-   Passed as a `qtl2` formatted genoprobs object. This option is for
    passing to `runQTL` and `mapQTL`. An example of this is shown below.

``` r
## Structure and example of qtl2 fomatted genoprobs
str(demo_probs)
#> List of 20
#>  $ 1 : num [1:10, 1:8, 1:4711] 1.48e-08 1.23e-08 1.04e-08 3.31e-08 2.94e-07 ...
#>   ..- attr(*, "dimnames")=List of 3
#>   .. ..$ : chr [1:10] "PB357.02_repA" "PB357.12_repA" "PB357.16_repA" "PB357.18_repA" ...
#>   .. ..$ : chr [1:8] "A" "B" "C" "D" ...
#>   .. ..$ : chr [1:4711] "1_3000000" "1_3041392" "1_3346528" "1_3651663" ...
#>  $ 2 : num [1:10, 1:8, 1:4709] 0.5 0.5 0.5 0.5 0.5 ...
#>   ..- attr(*, "dimnames")=List of 3
#>   .. ..$ : chr [1:10] "PB357.02_repA" "PB357.12_repA" "PB357.16_repA" "PB357.18_repA" ...
#>   .. ..$ : chr [1:8] "A" "B" "C" "D" ...
#>   .. ..$ : chr [1:4709] "2_3000000" "2_3038312" "2_3076624" "2_3510793" ...
#>  $ 3 : num [1:10, 1:8, 1:3811] 2.65e-09 5.42e-10 5.00e-01 1.32e-09 5.00e-01 ...
#>   ..- attr(*, "dimnames")=List of 3
#>   .. ..$ : chr [1:10] "PB357.02_repA" "PB357.12_repA" "PB357.16_repA" "PB357.18_repA" ...
#>   .. ..$ : chr [1:8] "A" "B" "C" "D" ...
#>   .. ..$ : chr [1:3811] "3_3000000" "3_3476607" "3_3953213" "3_4429820" ...
#>  $ 4 : num [1:10, 1:8, 1:3872] 3.51e-11 3.77e-11 3.51e-11 7.84e-11 6.94e-11 ...
#>   ..- attr(*, "dimnames")=List of 3
#>   .. ..$ : chr [1:10] "PB357.02_repA" "PB357.12_repA" "PB357.16_repA" "PB357.18_repA" ...
#>   .. ..$ : chr [1:8] "A" "B" "C" "D" ...
#>   .. ..$ : chr [1:3872] "4_3000000" "4_3039973" "4_3079946" "4_3119919" ...
#>  $ 5 : num [1:10, 1:8, 1:3837] 1.33e-18 7.70e-11 7.50e-20 2.88e-08 4.65e-09 ...
#>   ..- attr(*, "dimnames")=List of 3
#>   .. ..$ : chr [1:10] "PB357.02_repA" "PB357.12_repA" "PB357.16_repA" "PB357.18_repA" ...
#>   .. ..$ : chr [1:8] "A" "B" "C" "D" ...
#>   .. ..$ : chr [1:3837] "5_3000000" "5_3039177" "5_3333513" "5_3343035" ...
#>  $ 6 : num [1:10, 1:8, 1:3653] 4.30e-08 4.18e-14 4.30e-08 7.89e-18 7.89e-18 ...
#>   ..- attr(*, "dimnames")=List of 3
#>   .. ..$ : chr [1:10] "PB357.02_repA" "PB357.12_repA" "PB357.16_repA" "PB357.18_repA" ...
#>   .. ..$ : chr [1:8] "A" "B" "C" "D" ...
#>   .. ..$ : chr [1:3653] "6_3000000" "6_3040509" "6_3407400" "6_3774290" ...
#>  $ 7 : num [1:10, 1:8, 1:4006] 4.77e-18 5.27e-13 2.61e-15 1.00e-15 5.96e-18 ...
#>   ..- attr(*, "dimnames")=List of 3
#>   .. ..$ : chr [1:10] "PB357.02_repA" "PB357.12_repA" "PB357.16_repA" "PB357.18_repA" ...
#>   .. ..$ : chr [1:8] "A" "B" "C" "D" ...
#>   .. ..$ : chr [1:4006] "7_3000000" "7_3035905" "7_3071810" "7_3183355" ...
#>  $ 8 : num [1:10, 1:8, 1:3387] 9.17e-21 1.01e-09 2.96e-08 5.00e-01 3.00e-08 ...
#>   ..- attr(*, "dimnames")=List of 3
#>   .. ..$ : chr [1:10] "PB357.02_repA" "PB357.12_repA" "PB357.16_repA" "PB357.18_repA" ...
#>   .. ..$ : chr [1:8] "A" "B" "C" "D" ...
#>   .. ..$ : chr [1:3387] "8_3000000" "8_3037619" "8_3075238" "8_3402939" ...
#>  $ 9 : num [1:10, 1:8, 1:3414] 5.00e-01 1.24e-08 5.00e-01 1.41e-08 7.55e-08 ...
#>   ..- attr(*, "dimnames")=List of 3
#>   .. ..$ : chr [1:10] "PB357.02_repA" "PB357.12_repA" "PB357.16_repA" "PB357.18_repA" ...
#>   .. ..$ : chr [1:8] "A" "B" "C" "D" ...
#>   .. ..$ : chr [1:3414] "9_3000000" "9_3035823" "9_3071646" "9_3522909" ...
#>  $ 10: num [1:10, 1:8, 1:3450] 2.85e-12 5.00e-01 2.85e-12 2.85e-12 3.50e-11 ...
#>   ..- attr(*, "dimnames")=List of 3
#>   .. ..$ : chr [1:10] "PB357.02_repA" "PB357.12_repA" "PB357.16_repA" "PB357.18_repA" ...
#>   .. ..$ : chr [1:8] "A" "B" "C" "D" ...
#>   .. ..$ : chr [1:3450] "10_3000000" "10_3037150" "10_3074301" "10_3111451" ...
#>  $ 11: num [1:10, 1:8, 1:3796] 1.49e-10 1.65e-02 5.00e-01 4.84e-04 2.87e-04 ...
#>   ..- attr(*, "dimnames")=List of 3
#>   .. ..$ : chr [1:10] "PB357.02_repA" "PB357.12_repA" "PB357.16_repA" "PB357.18_repA" ...
#>   .. ..$ : chr [1:8] "A" "B" "C" "D" ...
#>   .. ..$ : chr [1:3796] "11_3000000" "11_3031502" "11_3063005" "11_3094507" ...
#>  $ 12: num [1:10, 1:8, 1:3124] 2.37e-01 1.42e-11 1.42e-11 1.18e-11 1.65e-09 ...
#>   ..- attr(*, "dimnames")=List of 3
#>   .. ..$ : chr [1:10] "PB357.02_repA" "PB357.12_repA" "PB357.16_repA" "PB357.18_repA" ...
#>   .. ..$ : chr [1:8] "A" "B" "C" "D" ...
#>   .. ..$ : chr [1:3124] "12_3000000" "12_3038001" "12_3076003" "12_3437198" ...
#>  $ 13: num [1:10, 1:8, 1:3229] 1.08e-10 5.00e-01 9.70e-11 5.00e-01 5.00e-01 ...
#>   ..- attr(*, "dimnames")=List of 3
#>   .. ..$ : chr [1:10] "PB357.02_repA" "PB357.12_repA" "PB357.16_repA" "PB357.18_repA" ...
#>   .. ..$ : chr [1:8] "A" "B" "C" "D" ...
#>   .. ..$ : chr [1:3229] "13_3000000" "13_3396518" "13_3793036" "13_3803732" ...
#>  $ 14: num [1:10, 1:8, 1:3019] 5.00e-01 5.00e-01 4.97e-01 1.34e-08 5.66e-06 ...
#>   ..- attr(*, "dimnames")=List of 3
#>   .. ..$ : chr [1:10] "PB357.02_repA" "PB357.12_repA" "PB357.16_repA" "PB357.18_repA" ...
#>   .. ..$ : chr [1:8] "A" "B" "C" "D" ...
#>   .. ..$ : chr [1:3019] "14_3000000" "14_3040840" "14_3081681" "14_3122521" ...
#>  $ 15: num [1:10, 1:8, 1:2761] 9.46e-11 6.94e-16 2.28e-13 5.00e-01 5.00e-01 ...
#>   ..- attr(*, "dimnames")=List of 3
#>   .. ..$ : chr [1:10] "PB357.02_repA" "PB357.12_repA" "PB357.16_repA" "PB357.18_repA" ...
#>   .. ..$ : chr [1:8] "A" "B" "C" "D" ...
#>   .. ..$ : chr [1:2761] "15_3000000" "15_3036880" "15_3073761" "15_3110641" ...
#>  $ 16: num [1:10, 1:8, 1:2688] 5.00e-01 5.00e-01 5.00e-01 1.76e-10 2.59e-10 ...
#>   ..- attr(*, "dimnames")=List of 3
#>   .. ..$ : chr [1:10] "PB357.02_repA" "PB357.12_repA" "PB357.16_repA" "PB357.18_repA" ...
#>   .. ..$ : chr [1:8] "A" "B" "C" "D" ...
#>   .. ..$ : chr [1:2688] "16_3000000" "16_3035541" "16_3071083" "16_3106624" ...
#>  $ 17: num [1:10, 1:8, 1:2873] 8.69e-21 2.36e-20 1.64e-22 5.00e-01 2.36e-20 ...
#>   ..- attr(*, "dimnames")=List of 3
#>   .. ..$ : chr [1:10] "PB357.02_repA" "PB357.12_repA" "PB357.16_repA" "PB357.18_repA" ...
#>   .. ..$ : chr [1:8] "A" "B" "C" "D" ...
#>   .. ..$ : chr [1:2873] "17_3000000" "17_3032196" "17_3064392" "17_3375528" ...
#>  $ 18: num [1:10, 1:8, 1:2588] 5.00e-01 4.86e-13 4.75e-13 5.00e-01 5.00e-01 ...
#>   ..- attr(*, "dimnames")=List of 3
#>   .. ..$ : chr [1:10] "PB357.02_repA" "PB357.12_repA" "PB357.16_repA" "PB357.18_repA" ...
#>   .. ..$ : chr [1:8] "A" "B" "C" "D" ...
#>   .. ..$ : chr [1:2588] "18_3000000" "18_3416510" "18_3833020" "18_4249530" ...
#>  $ 19: num [1:10, 1:8, 1:2434] 5.00e-01 5.84e-20 2.83e-01 3.73e-20 2.93e-01 ...
#>   ..- attr(*, "dimnames")=List of 3
#>   .. ..$ : chr [1:10] "PB357.02_repA" "PB357.12_repA" "PB357.16_repA" "PB357.18_repA" ...
#>   .. ..$ : chr [1:8] "A" "B" "C" "D" ...
#>   .. ..$ : chr [1:2434] "19_3000000" "19_3024080" "19_3048159" "19_3072239" ...
#>  $ X : num [1:10, 1:8, 1:3643] 5.13e-20 5.13e-20 9.40e-01 1.90e-02 2.77e-02 ...
#>   ..- attr(*, "dimnames")=List of 3
#>   .. ..$ : chr [1:10] "PB357.02_repA" "PB357.12_repA" "PB357.16_repA" "PB357.18_repA" ...
#>   .. ..$ : chr [1:8] "A" "B" "C" "D" ...
#>   .. ..$ : chr [1:3643] "X_3000000" "X_3046778" "X_3093557" "X_3140335" ...
#>  - attr(*, "is_x_chr")= logi [1:20] FALSE FALSE FALSE FALSE FALSE FALSE ...
#>  - attr(*, "crosstype")= chr "DO"
#>  - attr(*, "alleles")= chr [1:8] "A" "B" "C" "D" ...
#>  - attr(*, "alleleprobs")= logi TRUE
#>  - attr(*, "class")= chr [1:2] "calc_genoprob" "list"
```

``` r
demo_probs[[1]][1:2,,1:2]
#> , , 1_3000000
#> 
#>                          A   B            C            D            E
#> PB357.02_repA 1.477899e-08 0.5 1.217364e-08 1.370868e-08 1.199720e-08
#> PB357.12_repA 1.232415e-08 0.5 5.000000e-01 1.604246e-08 3.499095e-09
#>                          F            G            H
#> PB357.02_repA 5.000000e-01 4.331041e-09 1.313560e-08
#> PB357.12_repA 1.478095e-08 4.742920e-09 1.118048e-08
#> 
#> , , 1_3041392
#> 
#>                          A   B            C            D            E
#> PB357.02_repA 1.477899e-08 0.5 1.217364e-08 1.370868e-08 1.199720e-08
#> PB357.12_repA 1.232415e-08 0.5 5.000000e-01 1.604246e-08 3.499095e-09
#>                          F            G            H
#> PB357.02_repA 5.000000e-01 4.331041e-09 1.313560e-08
#> PB357.12_repA 1.478095e-08 4.742920e-09 1.118048e-08
```

### Quantified Data

Although the variable is named `expr_mats` any quantified molecular data
can be passed to this variable. They can be passed as a vector of paths
to the data (in the same order as tissues are being passed), or as a
named list of matrices (named with the tissue type associated with the
counts). These should be filtered, batch corrected and normalized data
with the phenotype in rows and the samples in columns. See the attached
`demo_counts` for an example of formatting.

``` r
## Information and partial example of quantified data
class(demo_counts)
#> [1] "matrix" "array"
```

``` r
head(demo_counts)
#>                    PB357.02_repA PB357.12_repA PB357.16_repA PB357.18_repA
#> ENSMUSG00000030402      90.27620      47.49881      65.92215      60.68472
#> ENSMUSG00000025473      55.33058      61.24847      59.68627      26.61611
#> ENSMUSG00000029246    2333.74400    5135.62454    2341.10901    2701.05069
#> ENSMUSG00000051413     218.41017     162.49593     209.34736     223.57530
#> ENSMUSG00000026356    3244.11907   11302.21694    2892.56671    4801.55973
#> ENSMUSG00000042505     317.42278     554.98610     316.24815     300.22969
#>                    PB357.19_repA PB357.20_repA PB357.21_repA PB357.22_repA
#> ENSMUSG00000030402       50.6948      43.98085      68.98128      57.27735
#> ENSMUSG00000025473       24.3725      45.02802      55.41119      37.55891
#> ENSMUSG00000029246     3481.3853    2611.54250    1373.65252    3451.91477
#> ENSMUSG00000051413      230.0764     158.12164     211.46720     217.84171
#> ENSMUSG00000026356     5616.3992    4653.59323    2312.56906    5754.02576
#> ENSMUSG00000042505      348.0393     347.65818     360.73816     294.83748
#>                    PB357.28_repA PB357.29_repA
#> ENSMUSG00000030402      42.11575      52.36949
#> ENSMUSG00000025473      14.03858      44.38093
#> ENSMUSG00000029246    3130.74325    2534.28661
#> ENSMUSG00000051413     161.91166     260.95984
#> ENSMUSG00000026356    6941.61141    4452.29448
#> ENSMUSG00000042505     440.81151     286.70078
```

### Marker Grid

This should be a grid of markers and the associated genomic locations in
bp and cM. The grid file can be passed as an object, or as a path to
where the file is located. There are two grid files attached to
`QTLretrievR`, a grid with ~75K pseudo-markers with positions updated
for mm38 and a grid with ~69K pseudo-markers with positions associated
with mm10. If your probabilities are derived from a MUGA/GIGAMUGA scan,
the markers can be interpolated onto either grid of your choice using
GBRS, or you can use a new marker grid to be used with your
probabilities. If making your own grid file, see `gridfile` to make sure
you have the correct format.

``` r
## example format for 75K marker grid
head(gridfile)
#>           chr     pos      cM      bp
#> 1_3000000   1 3000000 0.00001 3000000
#> 1_3039563   1 3039563 0.02001 3039563
#> 1_3079126   1 3079126 0.04001 3079126
#> 1_3118689   1 3118689 0.06001 3118689
#> 1_3158252   1 3158252 0.08001 3158252
#> 1_3197814   1 3197814 0.10001 3197814
```

### Metadata

At a minimum, the metadata file should contain a column for sample IDs
(`ID`) and columns for any covariates to be included in the peak and
effect calling. Any other pieces of information that may be relevant to
your analysis that you want to make sure is all in the same place should
also be included. A minimal example, `demo_meta`, is attached to
`QTLretrievR`.

``` r
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

### Annotations

Currently the annotations needs to contain at a minimum `gene`,
`symbol`, `chr`, `start`, and `end` columns. The `gene` column should
contain the phenotype information that appears in the row names of the
quantified data. In the attached example (`demo_annot`) the gene column
contains ensembl gene IDs, these could be traded out for ensembl protein
IDs or any other phenotype identifier.

``` r
## Example annotation file - Note this file is all the genes for this biomart build
head(demo_annot)
#>                 gene symbol chr     start       end
#> 1 ENSMUSG00000000001  Gnai3   3 108107280 108146146
#> 2 ENSMUSG00000000028  Cdc45  16  18780447  18811987
#> 3 ENSMUSG00000000031    H19   7 142575529 142578143
#> 4 ENSMUSG00000000037  Scml2   X 161117193 161258213
#> 5 ENSMUSG00000000056   Narf  11 121237253 121255856
#> 6 ENSMUSG00000000058   Cav2   6  17281185  17289115
```

## Running as Separate Steps

If, instead of using the wrapper to run all steps at once, you wish to
run each step individually (calculating/formatting the genotype
probabilities, mapping/peak calling, mediation, effects) this can be
done with relative ease.

### Genoprobs

In order to convert the genoprobs from the GBRS output to the `qtl2`
format we run a function called `genoprobably`. It pulls in the relevant
files from the provided location and splits samples into the tissues
they are derived from (the names that you pass to `tissue` should be
included in the file name in some form), then puts these into a 3
dimensional format based on the markers. Once the genotype probabilities
are in the 3D format `probs_3d_to_qtl2` is used to convert them to the
qtl2 format. The below example shows how to run genoprobably with the
75k gridfile.

``` r
## Unevaluated Code Chunk
probs <- genoprobably(outfile = "../../vignette/gbrs_interpolated_probs.rds", 
                      gbrsFileLoc = "<path/to/gbrs/results>",
                      tissues = c("tissue1","tissue2"), 
                      gridFile = gridfile, 
                      save = "ro")
```

### Mapping and Peak Calling

The mapping and peak calling steps are combined in our pipeline. We use
the genoprobs to calculate kinship, and create the genetic and pysical
maps. Then using the quantified data and sample metadata we rankZ
normalize the counts, calculate covariates and finally do the peak
calling using `qtl2::scan1`.

``` r
## Unevaluated Code Chunk
## Passing objects:
map_peaks <- mapQTL(outdir = "../../vignette/", 
                    peaks_out = "mm39_peaks.rds",
                    map_out = "mm39_map.rds", 
                    genoprobs = probs, 
                    samp_meta = demo_meta,
                    expr_mats = list("mESC" = demo_counts), 
                    covar_factors = c("sex"),
                    gridFile = gridfile, 
                    biomart = demo_annot, 
                    save = "ro")

## Passing file locations
map_peaks <- mapQTL(outdir = "../../vignette/", 
                    peaks_out = "mm39_peaks.rds",
                    map_out = "mm39_map.rds", 
                    genoprobs = "path/to/saved/probs.rds",
                    samp_meta = "path/to/sample/metadata",
                    expr_mats = c("path/to/counts1","path/to/coounts/2"), 
                    covar_factors = c("sex"),
                    gridFile = "path/to/custom/grid",
                    biomart = "path/to/annotations",
                    save = "ro")
```

### Mediation

For mediation we use the peaks, mapping information (which includes the
rankZ normalized counts), and annotations to identify targets which in
turn allows us to identify phenotypes that act as mediators for peaks
above a suggestive LOD. This is done using
`intermediate::mediation.scan`.

``` r
## Unevaluated Code Chunk
med_res <- run_mediate(peaks = map_peaks$peaks_list, 
                       mapping = map_peaks$maps_list,
                       outdir = "../../vignette", 
                       biomart = demo_annot,
                       suggLOD - 7, 
                       med_out = "mm39_mediation.rds", 
                       save = "ro")
```

### Effects

For mediation we use the peaks, mapping information (which includes the
rankZ normalized counts), and annotations to identify targets which in
turn allows us to identify phenotypes that act as mediators for peaks
above a suggestive LOD. This is done using
`intermediate::mediation.scan`.

``` r
## Unevaluated Code Chunk
effects <- qtl_effects(mapping = map_peaks$maps_list, 
                       peaks = map_peaks$peaks_list,
                       suggLOD = 8, 
                       outdir = "../../vignette", 
                       outfile = "mm39_effects.rds", 
                       save = "ro")
```
