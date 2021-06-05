methodname ='Mod',  Com_Methodname_abbr='MO'
% % % methodname='QsSprs'  
% % % methodname='HC-PIN-log'    
Com_Methodname= methodname  
% % % % % % % % % 
netname    ='PPI'  
DisSetName ='Mesh'
interval_n_partition = 20 
% % % % % % % % % % % % % % % % % 
fdatestrForfile = datestr(now,'yyyy.mmm.dd-HH.MM.SS');  
fdataset        = ['data',filesep,'demoDataSet&PPICOM_ModCM_delta=0.2.mat']
load(fdataset) 
Matrix_gene_dis00      =  AdjGfD ; 
n_disgenes_eachdisease = sum(Matrix_gene_dis00,1)'; 
COMtable_ct            = COMtable;  

% % % % % % % % % % %    
rng(1111)
[n_gene_all, n_disease] = size(Matrix_gene_dis00);           
% % nCVTimes  = 2,   n_fold    = 5 ,   CVtype = [num2str(n_fold),'-foldCV'], MinSizeDisGeSet = n_fold ;  
nCVTimes  = 50,   n_fold    = 5 ,   CVtype = [num2str(n_fold),'FCV'], MinSizeDisGeSet = n_fold,  

dis_IDset = find(n_disgenes_eachdisease>=MinSizeDisGeSet); 
% dis_IDset = dis_IDset(1:2)   %%%for test only %%%%%%%%%%%    
n_dis =  length( dis_IDset )  

n_disease_in_Table   = length( dis_IDset   ); 
nCV_list             = zeros( n_disease_in_Table, 1 );  	
matAUROC_nCVTimes    = cell(nCVTimes,1);
matAUPRC_nCVTimes    = cell(nCVTimes,1);
matRec5_nCVTimes     = cell(nCVTimes,1);
matPrec5_nCVTimes    = cell(nCVTimes,1);
methodset_nCVTimes   = cell(nCVTimes,1);
parfor  i_cv = 1:nCVTimes    
    disp(['i_cv-',num2str(i_cv) ]) 
    %
    matAUROC    = [];
    matAUPRC    = [];
    matRec5     = []; 
    matPrec5    = []; 
    idx_res     = 0;     
    for ii_dis = 1:n_disease_in_Table
        tic
        Matrix_gene_dis_copy = Matrix_gene_dis00 ; 
        COM_Dataset          = COMtable_ct;  
        ID_dis               = dis_IDset(ii_dis);  
        % disp(['i_cv-',num2str(i_cv),'; ii_dis-',num2str(ii_dis),'; ID_dis-',num2str(ID_dis)]) 
        ac_gene_dis00 = Matrix_gene_dis_copy(:,ID_dis ); 
        idx_pos       = find( ac_gene_dis00 );  n_pos = length( idx_pos); 
        idx_neg       = find( ~ac_gene_dis00 ); n_neg = length( idx_neg); 
        n_fold_real   = min(n_fold, n_pos) ;   
        ind_fold_pos  = crossvalind('Kfold', n_pos, n_fold_real ) ; 
        ind_fold_neg  = crossvalind('Kfold', n_neg, n_fold_real ) ; 
        for i_fold = 1:n_fold_real 
            t1 = toc ; 
            % idx_pos_train = idx_pos(ind_fold_pos~=i_fold);
            idx_pos_test    = idx_pos(ind_fold_pos==i_fold);  n_pos_test =length(idx_pos_test); 
            %   
            idx_neg_test_WG        = idx_neg ;  n_neg_test_all = length(  idx_neg_test_WG  ) ;               
            idx_test_pos_neg_WG    = [idx_neg_test_WG; idx_pos_test ] ;  
            AdjGfD                       = Matrix_gene_dis_copy; 
            AdjGfD(idx_pos_test,ID_dis ) = 0 ;    

            %%%%%%%%%%%%% 
            TableScores1 = table; TableScores2= table; 
            % use only PPI 
            % [TableScores1, COM_Dataset ] = A_HyMM(COM_Dataset, AdjGfG,AdjGfD,[],     ID_dis );   
            % TableScores1.Properties.VariableNames = strcat(TableScores.Properties.VariableNames,'_HO') ; 

            % use disease similarity       
            [TableScores2, COM_Dataset ] = A_HyMM(COM_Dataset, AdjGfG,AdjGfD,AdjDfD, ID_dis );          
            % TableScores2.Properties.VariableNames = strcat(TableScores.Properties.VariableNames,'_HE') ; 

            %%%%%%%%%%%%%                              
            TableScores = [TableScores1,TableScores2]; 
            methodset = TableScores.Properties.VariableNames ;  
            %
            test_real = ac_gene_dis00(idx_test_pos_neg_WG);  
            [AUROCset, AUPRCset , Rec5set  , Prec5set ,n_method, methodnames] = getPerf(test_real,TableScores(idx_test_pos_neg_WG,:) ) ; 
            idx_res = idx_res +1 ; 
            matAUROC(idx_res,:) = AUROCset;
            matAUPRC(idx_res,:) = AUPRCset;
            matRec5(idx_res,:)  = Rec5set; 
            matPrec5(idx_res,:) = Prec5set; 
            methodset_nCVTimes{i_cv}= methodset;  
        end 
        %
        matAUROC_nCVTimes{i_cv} = mean(matAUROC,1);
        matAUPRC_nCVTimes{i_cv} = mean(matAUPRC,1);
        matRec5_nCVTimes{i_cv}  = mean(matRec5,1); 
        matPrec5_nCVTimes{i_cv} = mean(matPrec5,1); 
        % toc
    end 
    toc 
    disp('  ')
% 
end 
%  
matAUROC_nCVTimes = cat(1,matAUROC_nCVTimes{:});    
matAUPRC_nCVTimes = cat(1,matAUPRC_nCVTimes{:});    
matRec5_nCVTimes  = cat(1,matRec5_nCVTimes{:}); 
matPrec5_nCVTimes = cat(1,matPrec5_nCVTimes{:}); 
% 
matRESmean = [mean(matAUROC_nCVTimes,1);mean(matAUPRC_nCVTimes,1);mean(matRec5_nCVTimes,1);mean(matPrec5_nCVTimes,1) ]   
tbRESmean  = array2table(matRESmean, 'VariableNames', methodset_nCVTimes{1}, 'RowNames',{'AUROC','AUPRC','Rec5','Prec5'}) 
% 
%  save 
dir_results = 'results';  
if ~exist(dir_results,'dir'); mkdir(dir_results);end 
date_cmplt = datestr(now,'yyyy.mmm.dd-HH.MM.SS');
parastr    = sprintf('CVtype=%s_CVtime=%d_MSDGS%d', CVtype  ,  nCVTimes, MinSizeDisGeSet );   
outfile    = [dir_results,filesep,'ResPerf_HyMM_',Com_Methodname_abbr,'_Ninterval_',num2str(interval_n_partition),'_',netname,'_',DisSetName,'_',parastr,'_',date_cmplt,'.mat'] 
save([outfile],  'tbRESmean',   'date_cmplt'  , '-v7.3' )   ;   

     
     
    
