# SPAR_Paper_Figures_Code
This repository gives the full reproducible code for generating all Figures and Tables in Sparse Projected Averaged Regression to High-dimensional Data (see [Parzer, Vana-Guer and Filzmoser 2023](https://doi.org/10.48550/arXiv.2312.00130)).

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
