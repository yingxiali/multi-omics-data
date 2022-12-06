# multi-omics-data
This is the electronic appendix (R code) to the paper 
"Can combining many data types in multi-omics data lead to a worsening of predictive performance? A large-scale benchmark study".

#This repository contains 5 subfolders.

*The first one is "Data" subfolder.
The subfolder "Download" contains the R code 
whose function is to download the data needed for replicating the analyses presented in this paper from OpenML using the file "datset_ids.RData", which is also included in this repository.

*The second one is "JobScripts" subfolder.
It contains three R scripts "AnalysisCluster_1_4_5.R", "AnalysisCluster_2.R", "AnalysisCluster_3.R".
These allow for reproducing the benchmark study using the functions from the "Functions" subfolder.

"AnalysisCluster_1_4_5.R" performs the benchmark study for single block, the combinations of 4 blocks, and the combinations of 5 blocks.

"AnalysisCluster_2.R" performs the benchmark study for the combinations of 2 blocks.

"AnalysisCluster__3.R" performs the benchmark study for the combinations of 3 blocks.

*The third one  is "Functions" subfolder.
It also contains three R scripts whose labels contain "Functions_AnalysisCluster". These feature functions which are used for applying the various considered feature selection and classification methods.


*The fourth one is "Evaluation" subfolder.
It contains three R scripts"Evaluation_AnalysisCluster_fivemethods.R", "bootstrap analysis_ibrier.R" , and "bootstrap analysis_cindex.R".
"Evaluation_AnalysisCluster_fivemethods.R" evaluates the raw results to produce files, which are used in other scripts to produce the figures.
"bootstrap analysis" performs bootstrap analysis

*The fifth one is "Evaluation" subfolder.
 This subfolder allows to reproduce the figures shown in the paper and in the supplement.

Details on these approaches can be found in the paper.
