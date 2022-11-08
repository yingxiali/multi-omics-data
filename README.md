# multi-omics-data
This is the electronic appendix (R code) to the paper 
"Can combining many data types in multi-omics data lead to a worsening of predictive performance? A large-scale benchmark study".

#This repository contains 5 main types of files.

The first one is "down_data.R" whose function is to download the data needed for replicating the analyses presented in this paper from OpenML using the file "datset_ids.RData", which is also included in this repository.

The second one are the R files whose labels contain "AnalysisCluster". These allow for reproducing the benchmark study using the functions from the R files whose labels contain "Functions_AnalysisCluster".

The third one are the R files whose labels contain "Functions_AnalysisCluster". These feature functions which are used for applying the various considered feature selection and classification methods.

The fouth one are the R files whose labels contain "Evaluation_AnalysisCluster_fivemethods.R". It allows replication to evaluate our results.

The five one are the R files whose labels contain "figures.R". It allows to reproduce the graphs of our results.

#Further details on the different R scripts

"AnalysisCluster_1_4_5.R" performs the benchmark study for single block, the combinations of 4 blocks, and the combinations of 5 blocks.

"AnalysisCluster_2.R" performs the benchmark study for the combinations of 2 blocks.

"AnalysisCluster__3.R" performs the benchmark study for the combinations of 3 blocks.

Details on these approaches can be found in the paper.
