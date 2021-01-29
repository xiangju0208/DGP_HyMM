% % function  [ResScalartable, ResCurvetable ,ResTopKtable, n_gene,n_dis] = getResPerfFromScorelist(M_gene_dis_postestlabel,Score_gene_dis, topk_set, topk_AP)
function  [AUROCset, AURecallset , Recall50set  , n_gene,n_method, methodnames] = getResPerfFromScorelist(testsetlabel,TableScores )
    %% calculate AUC for one disease each time ....  % % % % % % % % % % % % % % % % % % % % % % %
    % input only list for testset, excluding traindata  
    methodnames = TableScores.Properties.VariableNames ; 
    [n_gene,n_method ] = size(TableScores) ;   %   
    AUROCset    = 0.5.*ones( 1, n_method    ); 
    AURecallset = 0.0.*ones( 1, n_method    );
    Recall50set = 0.0.*ones( 1, n_method    ); 
    top_k_vec = (1:200)' ; len_topk  = length( top_k_vec  ); 
     
    labels = testsetlabel;
    for ii=1:n_method  
        scores = double(  TableScores{:, ii }    );  
        N_pos  = nnz(  labels   );  
        % AUROC 
        [X,Y,~,AUC ] = perfcurve( labels , scores ,1 );          
        AUROCset(ii) = AUC ;    
          
        % top k 
        % top k    reacall   precision 
        top_k_vec_t = top_k_vec; top_k_vec_t(top_k_vec_t>n_gene)=n_gene;
        [X00,Y00 ] = perfcurve( labels , scores , 1,  'XCrit','tp+fp', 'YCrit','tp' ) ;   %%[X00,Y00 ]
        [~,Y] = getUniformXY(X00,Y00,top_k_vec_t, [],   [], 'linearperturbation' , 0 ) ;  
        topk=100; 
        AURecallset(ii)= sum(Y(top_k_vec_t<=topk)/N_pos) ;  
        Recall50set(ii)=  Y(top_k_vec_t==topk)/N_pos ;   
% %         
        
% %         topk=100;  
% %         [X00,Y00 ] = perfcurve( labels , scores , 1,  'XCrit','tp+fp', 'YCrit','tp' , 'XVals',[0:1:n_gene]) ;   
% %         AURecallset(ii)= sum(Y00(X00<=topk)/N_pos) ;   
% %         [X00(1:100)   Y00(1:100) ]
% %         plot(X00,Y00/N_pos)
% %         size((Y00(X00==50)/N_pos))
% %         Recall50set(ii)=  (Y00(X00==50)/N_pos) ;   
         
    end 
    
  

end
 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [XVals,YVals ] = getUniformXY(X,Y,XVals, funiform, y0, UniqueMethod , IsExtrap )  
% to facilitate the average curve, the x-axis is standardized 
% X,Y have the same length
% XVals specifies the values of x for output 
% Y_u: after taking the unique values X_u of X, average the Y values of the same x coordinate 
% funiform = @mean,   @max,    @min, or   @median, specifying how to obtain Y_u  
% YVals: according to the input XVals, obtained by linear interpolation  
% y0: if there does esits x=0, insert (x0,y0)
% this function is used to uniform the performance curves. 
% in many cases, UniqueMethod 'linearperturbation' without 'Extrap' is more suitable becase they are step-like. 
% in X,Y, the order of data is important information for the uniform
% operations, when 'linearperturbation' is used. It should be the increasing ordering of X and/or Y, similar to ROC
% and PRC data.  
if ~exist('XVals','var') || isempty( XVals  )  
    XVals = []; 
end

if ~exist('funiform','var') || isempty( funiform  )
    funiform = @mean; 
end

if ~exist('y0','var') || isempty( y0  )
    y0 = [] ; 
end

if ~exist('UniqueMethod','var') || isempty( UniqueMethod  )
    UniqueMethod = 'linearperturbation' ;  
end

if ~exist('IsExtrap','var') || isempty( IsExtrap  )
    IsExtrap = false ;  
end

X= X(:);
Y= Y(:);  

if nargin<3 || strcmp( XVals, 'unique') || isempty( XVals )
    [X_u,~,ic] = unique(X) ; 
    Y_u = zeros( size(X_u) ); 
    for ii = 1: length( X_u )
        % Y_u( ii ) = mean( Y(ic==ii) ) ;
        Y_u( ii ) = funiform( Y(ic==ii) ) ;
    end
    XVals = X_u;
    YVals = Y_u;   
    if nargout>2
        warning('Output1-2 are used for unique X; Output3-4 are unneccesary!');
    end
    return; 

elseif isnumeric(XVals) 
    XVals =XVals(:) ;
    switch UniqueMethod  
        case 'linearperturbation'  
            if  ~( any( X==0 ) ) && ~isempty( y0  )  % set y0 at x=0 if no point of x=0 
                X = [0; X ];
                Y = [y0;Y ];
            end
            % check if not extrap
            if ~IsExtrap
                X_max = max( X ) ;  
                X_min = min( X ) ;    
                if ~( all( X_min<=XVals )  && all( XVals <=X_max ) )
                    error('IsExtrap is false. The data cannot be extraped.');
                end            
            end
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %   
            X_rescale = max( abs(X) ); 
            if X_rescale~=0 
               if X_rescale ~=1     % normalization for more effective op 
                    X     = X./X_rescale;     
                    XVals = XVals./X_rescale;     
               end             
            else
                error('X_rescale: x_max==0'); 
            end
            % %   % % % % % % % % % % % % % % % % % % % % % % % %    
            pEPS     = 10^-12 ;
            length_X = length(X); 
            % % %             if length_X*pEPS>1/length(XVals)/10  
            if length_X*pEPS>(max(XVals)-min(XVals))/length(XVals)/10   
                warning(['List of X is too long: length(X)=',num2str(length_X)]); 
            elseif length_X<2 
                error(['List of X is too short: length(X)=',num2str(length_X)]); 
            end          
            %  
            X = X + [0:length_X-1]'.*pEPS;  % %  
            %          
            try 
                if IsExtrap
                    YVals = interp1(X ,Y ,XVals,'linear','extrap') ; 
                else
                    YVals = interp1(X ,Y ,XVals,'linear' ) ;             
                end 
            catch 
                disp( min(X) )
                disp( max(X) )
                disp( min(XVals) )
                disp( max(XVals) )
                [ min(XVals)<min(X)        max(X)<max(XVals)]
            end
            %
            if X_rescale~=1 
                XVals = XVals*X_rescale; 
            end 
            return; 
 
      
        otherwise
            error( 'UniqueMethod is wrong. It should be linearperturbation or unique.' );
    end    
    
else 
	error('Input parameter: XVals is wrong for getUniformXY(X,Y,XVals)!') ;     
end
 
end







