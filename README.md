# SPAR_Paper_Figures_Code
This repository gives the full reproducible code for generating all Figures and Tables in Sparse Projected Averaged Regression to High-dimensional Data (see [Parzer, Filzmoser and Vana-Guer 2024](https://doi.org/10.48550/arXiv.2312.00130)).
Published version: Sparse data-driven random projection in regression for high-dimensional data, [Parzer, Filzmoser and Vana-Guer 2025](https://doi.org/10.52933/jdssv.v5i5.138))

This repository consists of the following folders with described contents.

- Chaper2Experiments: R-scripts and .rds-files for the experiments and Figures in chapter/section 2 'Methods'
- data: .txt files for the rat data and a .mat file for the face angel data
- data_application: R-script for data applications applying all methods on all data sets multiple times and saving the resulting .rds file to the folder 'saved_results'; and another R-script reproducing the preprocessing of the face data in 'Compressed Gaussian Process for Manifold Regression' (Guhaniyogi and Dunson 2016) for performance comparison
- functions: 4 R-scripts, data_generation.R for defining a function generating data from a HD linear model, methods.R defining consistent wrapper functions for all considered methods, multi_assign.R to define an operator assigning multiple variables at once (by Daniel Kapla, TU Wien) and RPM_generation.R to define three different functions generating certain random projection matrices (Sparse, Sparse_CW, our proposed adapted Sparse_CW)
- generate_plots: R-scipts reading in .rds files from 'saved_results' and generating the plots and tables for the simulation study and the data applications and saving the plots as pdfs in 'plots'
- plots: all pdf Figures
- saved_results: .rds files produced from 'simulation' or 'data_application' folders
- simulations: R-script for simulation study applying all methods on all simulation settings multiple times and saving the resulting .rds file to the folder 'saved_results'
- TARP-master: R-code from 'Targeted Random Projection for Prediction From High-Dimensional Features' [(Mukhopadhyay and Dunson 2020)](https://github.com/david-dunson/TARP) adapted to return the estimated beta regression coefficient

To enhance reproducibility, the following output was obtained from `sessionInfo()` after executing the header of the file 'simulations_SPAR_CV_022024.R'.

R version 4.2.1 (2022-06-23)

Platform: aarch64-apple-darwin20 (64-bit)

Running under: macOS 14.3.1

Matrix products: default

LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib

locale:

[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:

[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:

 [1] spar_4.0.0        stringr_1.5.1     robustHD_0.8.0    robustbase_0.99-0 perry_0.3.1       ggplot2_3.5.1     SplitReg_1.0.2    MASS_7.3-60  
 
 [9] SIS_0.8-8         glmnet_4.1-8      pls_2.8-1         Matrix_1.6-1.1    ROCR_1.0-11       dplyr_1.1.4       tidyr_1.3.1       foreach_1.5.2    

loaded via a namespace (and not attached):

 [1] Rcpp_1.0.13-1    DEoptimR_1.0-11  pillar_1.9.0     compiler_4.2.1   iterators_1.0.14 tools_4.2.1      lifecycle_1.0.4  tibble_3.2.1     gtable_0.3.6   
 
[10] lattice_0.20-45  pkgconfig_2.0.3  rlang_1.1.4      cli_3.6.3        rstudioapi_0.13  withr_3.0.2      generics_0.1.3   vctrs_0.6.5      grid_4.2.1   

[19] tidyselect_1.2.1 glue_1.8.0       R6_2.5.1         fansi_1.0.6      Rdpack_2.6.2     survival_3.3-1   pacman_0.5.1     purrr_1.0.2      magrittr_2.0.3  

[28] rbibutils_2.3    scales_1.3.0     codetools_0.2-18 splines_4.2.1    colorspace_2.1-1 shape_1.4.6.1    ncvreg_3.14.1    utf8_1.2.4       stringi_1.8.4 

[37] munsell_0.5.1
