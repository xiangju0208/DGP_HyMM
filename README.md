# HyMM
HyMM: Hybrid method for disease-gene prediction by integrating multiscale module structures


## Requirements
Matlab 2016 or above   


## Codes 
#main_HyMM.m: cross-validation code.  <br>
This code allows parallel execution. You can change "parfor" to "for" to cancel parallel execution  <br>
 <br>
#A_DGP_HyMM_ByRank.m: the recommended HyMM algorithms in the study. <br>
A_DGP_HyMM(COM_Dataset, AdjGfG,AdjGfD,AdjDfD, DisIDset, plus_method_set, RankMergeMethod  )   
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
Xiang, et el. HyMM: Hybrid method for disease-gene prediction by integrating multiscale module structures


## contact<br>
Email: xiang.ju@foxmail.com 


 
