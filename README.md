# SPAR_Paper_Figures_Code
This repository gives the full reproducible code for generating all Figures and Tables in Sparse Projected Averaged Regression to High-dimensional Data (see [Parzer, Filzmoser and Vana-Guer 2024](https://doi.org/10.48550/arXiv.2312.00130)).

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

To enhance reproducibility, the following output was obtained from `sessionInfo()` after executing the header of the file 'simulations_SPAR_CV_022025.R'.

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

 [1] R.utils_2.12.2       R.oo_1.25.0          R.methodsS3_1.8.2    SLOPE_0.5.2          kableExtra_1.3.4     knitr_1.44          

 [7] ggrepel_0.9.4        owl_0.1.1            randomForest_4.7-1.1 SPAR_1.1.1           stringr_1.5.1        robustHD_0.8.0      

[13] robustbase_0.99-0    perry_0.3.1          ggplot2_3.5.1        SplitReg_1.0.2       MASS_7.3-60          SIS_0.8-8           

[19] glmnet_4.1-8         pls_2.8-1            Matrix_1.6-1.1       ROCR_1.0-11          dplyr_1.1.4          tidyr_1.3.1         

[25] foreach_1.5.2       

loaded via a namespace (and not attached):

 [1] nlme_3.1-163        xts_0.13.2          doParallel_1.0.17   webshot_0.5.3       RColorBrewer_1.1-3  httr_1.4.7         

 [7] tools_4.2.1         utf8_1.2.4          R6_2.6.1            colorspace_2.1-1    nnet_7.3-17         withr_3.0.2        

[13] tidyselect_1.2.1    greybox_2.0.2       curl_5.2.0          compiler_4.2.1      textshaping_0.3.6   cli_3.6.4          

[19] rvest_1.0.2         pacman_0.5.1        xml2_1.3.3          labeling_0.4.3      tsutils_0.9.4       tseries_0.10-55    

[25] scales_1.3.0        lmtest_0.9-40       DEoptimR_1.0-11     fracdiff_1.5-3      quadprog_1.5-8      systemfonts_1.0.4  

[31] digest_0.6.37       rmarkdown_2.25      svglite_2.1.0       pkgconfig_2.0.3     htmltools_0.5.7     fastmap_1.1.1      

[37] highr_0.9           rlang_1.1.5         TTR_0.24.4          rstudioapi_0.13     quantmod_0.4.26     ncvreg_3.14.1      

[43] shape_1.4.6.1       ggh4x_0.2.8         farver_2.1.2        generics_0.1.3      zoo_1.8-12          magrittr_2.0.3     

[49] texreg_1.39.4       Rcpp_1.0.13-1       munsell_0.5.1       lifecycle_1.0.4     stringi_1.8.4       forecast_8.23.0    

[55] grid_4.2.1          lattice_0.20-45     splines_4.2.1       pillar_1.10.1       codetools_0.2-18    pkgload_1.4.0      

[61] urca_1.3-4          glue_1.8.0          evaluate_1.0.0      BiocManager_1.30.22 vctrs_0.6.5         nloptr_2.0.3       

[67] gtable_0.3.6        purrr_1.0.2         xfun_0.42           xtable_1.8-4        MAPA_2.0.7          pracma_2.3.8       

[73] ragg_1.2.5          survival_3.3-1      viridisLite_0.4.2   timeDate_4021.104   tibble_3.2.1        iterators_1.0.14   

[79] statmod_1.5.0       smooth_4.0.2 
