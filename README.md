#  Code to the benchmark study based on multi-omics data by Li et al.
This is the electronic appendix (R code) to the article 
"Can combining many data types in multi-omics data lead to a worsening of predictive performance? A large-scale benchmark study" by
Yingxia Li, Ulrich Mansmann, and Roman Hornung.


#This repository contains 5 subfolders.


*The first one is the "Data" subfolder.
The function of the R script "down_data.R" is to download the data needed for replicating the analyses presented in this paper from OpenML using the file "datset_ids.RData", which is also included in this subfolder.


*The second one is the "JobScripts" subfolder.
It contains three R scripts "AnalysisCluster_1_4_5.R", "AnalysisCluster_2.R", and "AnalysisCluster_3.R".
These allow for reproducing the benchmark study using the functions from the "Functions" subfolder.

"AnalysisCluster_1_4_5.R" performs the benchmark study for single blocks, combinations of 4 blocks, and combinations of 5 blocks.

"AnalysisCluster_2.R" performs the benchmark study for combinations of 2 blocks.

"AnalysisCluster_3.R" performs the benchmark study for combinations of 3 blocks.


*The third one is the "Functions" subfolder.
It also contains three R scripts whose labels contain "Functions_AnalysisCluster". These feature functions which are used for applying the different prediction methods to the different combinations.


*The fourth one is the "Evaluations" subfolder.
It contains six R scripts "Evaluation_AnalysisCluster_fivemethods.R", "bootstrap analysis_ibrier.R", "bootstrap analysis_cindex.R", figures_1_S1_S2_S3.R, figures_2_S4.R, and figures_S5_S6_1.R.
"Evaluation_AnalysisCluster_fivemethods.R" evaluates the raw results to produce files, which are used in other scripts to produce the figures.
The R scripts whose labels contain "bootstrap_analysis" perform the bootstrap analysis.
The R scripts whose labels contain "figures" allow to reproduce the figures shown in the paper and in the supplement.


*The fifth one is the "Results" subfolder.
In the "rda_files" subfolder, the files "scenariogrid1.Rda", "scenariogrid2.Rda", and "scenariogrid3.Rda",  were generated by the R scripts contained in the folder "JobScripts"; for details, see these R scripts;
 The files "resultsumsum.RData", "resultsum.RData", CI_cindex.xlsx, and CI_ibrier.xlsx  were generated by the R scripts contained in the folder "Evaluations"; for details, see these R scripts.
 In the "figures" subfolder, it contains all figures in the paper and in the supplement and were produced by the R scripts contained in the folder "Evaluations"; for details, see these R scripts.


Details on the benchmark study can be found in the paper.
