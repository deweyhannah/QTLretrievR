# Package index

### Core Functions

Primary functions for peak calling, haplotype effect analysis, and
mediation

- [`genoprobably()`](https://deweyhannah.github.io/QTLretrievR/reference/genoprobably.md)
  : Convert GBRS tsv genome probabilities to genoprobs format.
- [`mugaprobs()`](https://deweyhannah.github.io/QTLretrievR/reference/mugaprobs.md)
  : Convert and process MUGA probabilities to qtl2 format
- [`runQTL()`](https://deweyhannah.github.io/QTLretrievR/reference/runQTL.md)
  : Wrapper function to generate mapping, peaks, mediation, and effects
  data
- [`mapQTL()`](https://deweyhannah.github.io/QTLretrievR/reference/mapQTL.md)
  : Generate mapping data and peaks for QTL analysis
- [`gxeQTL()`](https://deweyhannah.github.io/QTLretrievR/reference/gxeQTL.md)
  : Generate mapping data and peaks for G x E QTL analysis
- [`LOD_thld()`](https://deweyhannah.github.io/QTLretrievR/reference/LOD_thld.md)
  : Determine significant and suggestive LOD thresholds
- [`qtl_effects()`](https://deweyhannah.github.io/QTLretrievR/reference/qtl_effects.md)
  : Calculate founder effects for significant peaks
- [`modiFinder()`](https://deweyhannah.github.io/QTLretrievR/reference/modiFinder.md)
  : Prepare and run mediation for a set of QTL peaks
- [`multi_modiFinder()`](https://deweyhannah.github.io/QTLretrievR/reference/multi_modiFinder.md)
  : Prepare and run between phenotype mediation for a set of QTL peaks
  with provided expression data.

### Plotting Functions

Functions used for plotting and visualizations

- [`plot_eqtlmap()`](https://deweyhannah.github.io/QTLretrievR/reference/plot_eqtlmap.md)
  : Plot QTL maps (peak vs gene)
- [`transbands()`](https://deweyhannah.github.io/QTLretrievR/reference/transbands.md)
  : Identify distal hotspots above a suggestive and significant LOD
  score.
- [`hsHapEffects()`](https://deweyhannah.github.io/QTLretrievR/reference/hsHapEffects.md)
  : Map founder haplotype effects for a given hotspot.
- [`medPlot_hotSpot()`](https://deweyhannah.github.io/QTLretrievR/reference/medPlot_hotSpot.md)
  : Identify top mediators within a hotspot (molecular QTL)
- [`medPlot_single()`](https://deweyhannah.github.io/QTLretrievR/reference/medPlot_single.md)
  : Identify mediators for provided features +/- x Mb from a given
  position on a given chromosome
- [`hsPeakPlot()`](https://deweyhannah.github.io/QTLretrievR/reference/hsPeakPlot.md)
  : Plot an individual peak or Principal Component of a hotspot.
- [`tailWag()`](https://deweyhannah.github.io/QTLretrievR/reference/tailWag.md)
  : Create a plot with founder effects for a given phenotype
- [`peak_plot()`](https://deweyhannah.github.io/QTLretrievR/reference/peak_plot.md)
  : Plot peaks associated with specific genes.

### Attached Data

Datasets built into QTLretrievR for examples or use in mapping

- [`demo_counts`](https://deweyhannah.github.io/QTLretrievR/reference/demo_counts.md)
  : Demo Counts
- [`demo_annot`](https://deweyhannah.github.io/QTLretrievR/reference/demo_meta.md)
  : Demo Metadata
- [`demo_probs`](https://deweyhannah.github.io/QTLretrievR/reference/demo_probs.md)
  : Demo Probabilities
- [`demo_annot`](https://deweyhannah.github.io/QTLretrievR/reference/demo_annot.md)
  : Demo Annotations
- [`annot_105`](https://deweyhannah.github.io/QTLretrievR/reference/annot_105.md)
  : Annotations file
- [`gridfile`](https://deweyhannah.github.io/QTLretrievR/reference/gridfile.md)
  : Grid Map
- [`gridfile69k`](https://deweyhannah.github.io/QTLretrievR/reference/gridfile69k.md)
  : Grid Map
- [`gridfile_mini`](https://deweyhannah.github.io/QTLretrievR/reference/gridfile_mini.md)
  : Grid Map
- [`gridfile_GM`](https://deweyhannah.github.io/QTLretrievR/reference/gridfile_GM.md)
  : Grid Map
- [`gridfile_MM`](https://deweyhannah.github.io/QTLretrievR/reference/gridfile_MM.md)
  : Grid Map
- [`gridfile_quilt`](https://deweyhannah.github.io/QTLretrievR/reference/gridfile_quilt.md)
  : Grid Map

### Helper Functions

Functions run during mapQTL that may be useful after the factt

- [`annotatePeaks()`](https://deweyhannah.github.io/QTLretrievR/reference/annotatePeaks.md)
  : Interpolate peak positions and attach annotations
- [`interp_bp()`](https://deweyhannah.github.io/QTLretrievR/reference/interp_bp.md)
  : interp_bp
- [`interp_cM()`](https://deweyhannah.github.io/QTLretrievR/reference/interp_cM.md)
  : interp_cM
- [`check_data()`](https://deweyhannah.github.io/QTLretrievR/reference/check_data.md)
  : Checking the inputted data for QTLretrievR
- [`split_map()`](https://deweyhannah.github.io/QTLretrievR/reference/split_map.md)
  : Background functions for QTLretrievR that will be called by the
  primary functions
