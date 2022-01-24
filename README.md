# HyMM: Hybrid method for disease-gene prediction by integrating multiscale module structures
Motivation: Identifying disease-related genes is an important issue in computational biology. Module structure widely exists in biomolecule networks, and complex diseases are usually thought to be caused by perturbations of local neighborhoods in the networks, which can provide useful insights for the study of disease-related genes. However, the mining and effective utilization of the module structure are still challenging in such issues as disease-gene prediction.    
Results: We propose a hybrid disease-gene-prediction method integrating multiscale module structure (HyMM), which can utilize multiscale information from local to global structure to more effectively predict disease-related genes. HyMM ex-tracts module partitions from local to global scales by multiscale modularity optimization with exponential sampling, and estimates the disease relatedness of genes in partitions by the abundance of disease-related genes within modules. Then, a probabilistic model for integration of gene rankings is designed in order to integrate multiple predictions derived from mul-tiscale module partitions and network propagation, and a parameter estimation strategy based on functional information is proposed to further enhance HyMM’s predictive power. By a series of experiments, we reveal the importance of module partitions at different scales, and verify the stable and good performance of HyMM compared with eight other state-of-the-arts and its further performance improvement derived from the parameter estimation. 
Conclusions: The results confirm that HyMM is an effective framework for integrating multiscale module structure to en-hance the ability to predict disease-related genes, which may provide useful insights for the study of multiscale module structure and its application in such issues as disease-gene prediction.  


## Requirements
Matlab 2016 or above   


## Codes 
#main_HyMM.m: cross-validation code.  <br>
This code allows parallel execution. You can change "parfor" to "for" to cancel parallel execution  <br>
 <br>
#A_HyMM.m: the recommended HyMM algorithm in the study. <br>
A_HyMM(COM_Dataset, AdjGfG,AdjGfD,AdjDfD, DisIDset, plus_method_set, RankMergeMethod  )   
% Input:  <br>
% COM_Dataset is a table recording the partition matrices; <br>
% AdjGfG: associatins between Genes (G) and Genes (G)   <br>
% AdjGfD: associatins between Genes (G) and Diseases (D)  <br>
% AdjDfD: associatins between Diseases (D) and Diseases (G)  <br>
% DisIDset: disease index  <br>
% plus_method_set: baseline algorithms. If plus_method_set is given, the results of baseline algorithms will be output.   <br>
% e.g. plus_method_set = {'RWRH'};  <br>
% RankMergeMethod: aggregation method <br>
% Ouput: <br>
% TableScores: a table whos variable record the scores of genes. <br>
% COM_Dataset: multiscale partitions that are preprocessed, facilating the usage of partition information in the latter, e.g. for cross-validation.  <br>


## Dataset
A dataset is located in the directory: data/demoDataSet&PPICOM_ModCM_delta=0.2.mat<br>
This dataset includes: <br>
(1) disease-gene associations, disease-disease associations and gene-gene associations;  <br>
(2) multiscale module partition matrices. <br>


## Results 
The results will be automatically saved into the directory: results.  

## cite
If you use HyMM in your research, please cite: <br>
Xiang, et al., HyMM: Hybrid method for disease-gene prediction by integrating multiscale module structure, bioRxiv (2021), doi: 10.1101/2021.04.30.442111.


## contact<br>
Email: xiang.ju@foxmail.com 


 
