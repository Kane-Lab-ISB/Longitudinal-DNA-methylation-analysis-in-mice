# Longitudinal-DNA-methylation-analysis-in-mice
Raw DNAm sequencing data are deposited at in the Gene Expression Omnibus (GEO) under the accession number GSE316651

Process data, discovery_cohort_DNAm_data.RData (contains methylation matrix data for the discovery cohort) and validation_cohort_DNAm_data.RData (contains methylation matrix data for the validation cohort), can be found at 10.5281/zenodo.21692964.

methy_prob_loci.RData contains probe loci matrix, and DNA_methylation_functions.R contains functions used in the analysis

R version 4.4.2 (2024-10-31)
Platform: aarch64-apple-darwin20
Running under: macOS 26.2

other attached packages:
 [1] doParallel_1.0.17        iterators_1.0.14         foreach_1.5.2            randomForest_4.7-1.2     survminer_0.5.0         
 [6] survival_3.8-3           reshape2_1.4.4           viridis_0.6.5            viridisLite_0.4.2        lmerTest_3.1-3          
[11] lme4_1.1-36              Matrix_1.7-1             gtsummary_2.1.0          gt_0.11.1                ComplexUpset_1.3.3      
[16] caret_7.0-1              lattice_0.22-6           lubridate_1.9.4          forcats_1.0.0            stringr_1.5.1           
[21] purrr_1.0.2              readr_2.1.5              tibble_3.2.1             tidyverse_2.0.0          readxl_1.4.3            
[26] tidyr_1.3.1              pheatmap_1.0.12          ggpubr_0.6.0             ggfortify_0.4.17         dplyr_1.1.4             
[31] variancePartition_1.36.3 BiocParallel_1.40.0      limma_3.62.2             ggplot2_3.5.1           
