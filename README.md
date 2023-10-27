# multi_season_analysis
Multi-season analysis of Field Scanalyzer, soil water content, and AZMET information.

Each file does the following:
- azmet.ipynb: Visualizes AZMET data across multiple seasons. To download AZMET data, [click here](https://ag.arizona.edu/azmet/06.htm).
- linear_mixed_modeling.R: Runs linear mixed models (LMMs) and extracts best linear unbiased predictors (BLUPs) and heritability. **Please note this has to be run before running the scripts below.**
- blups_by_genotype.ipynb: Plots BLUPs per genotype.
- heritability.ipynb: Plots heritability values
- trait_fpca.ipynb: Run functional principal component analysis (FPCA) on phenotype information.
