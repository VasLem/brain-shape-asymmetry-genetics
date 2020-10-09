function [P_sec,T_sec,options] = NPC(T,comb_funct,options)
%Nonparametric Combination Methodology.
%P_sec = NPC(T,comb_funct) nonparametrically combines hypotheses
%producing a multivariate test. T is a (B+1)-by-m-by-k matrix of the
%permutation space of the statistics or p-values where B is the number of
%permutations, M the number of variables and K the number of strata. Each
%row corresponds to p-values or statistics related to a random permutation
%for the m variables. The last one is related to the observed data. This
%matrix can be obtained from functions such as NP_1s, NP_2s, NP_cs or similar.
%The default combining function is the Tippett's.
%P_sec is the (B+1) length P-values matrix with one dimension less than T. The
%P-Values refer to the combination of hypotheses along the dimension
%of T specified by OPTIONS.Combdims (see below). By default, the dimension
%to combine is the greatest dimension of T. The last row of P is related to
%the observed sample.
%
%P_sec = NPC(T,comb_funct) comb_funct specify the combination function to
%use. Possible settings for comb_funct are:
%
%           'Direct'      (or 'D')  combining function of t-statistic,
%           'Max-t'       (or 'M')  combining function of t-statistic,
%           'Fisher'      (or 'F')  combining function of p-values,
%           'Liptak'      (or 'L')  combining function of p-values,
%           'LiptakLog'   (or 'LL') combining function of p-values,
%           'Mahalanobis' (or 'Ma') combining function of p-values 
%                                   (only for univariate bidirectional
%                                   alternatives, but T computed using
%                                   unidirectional alternatives)
%           'Tippett'     (or 'T')  combining function of p-values (min-p).
%
%
%If comb_funct is a cell array, it performs an ITERATED PROCEDURE, for
%example:
%          
%       comb_funct{1}={'Tippett' [1 2]; 'Fisher' [1 2]; 'Liptak' [1 2]}
%       comb_funct{2}={'Fisher' [1 1 1]}
%   
%performs a two-level iterated procedure: at the first level (cell), it
%combines 2 partial p-values using Tippett's, Fisher's and Liptak's
%combining function. It gives double weight to the second hypothesis.
%The second level (cell)makes use of Fishe's combining function. Use the
%second column of the cells to assign weights (equally weighted hypotheses
%by default). By default, Tippet's combining function is used.
%
%P_sec = NPC(...,'OPTIONS') Possible settings for options are:
%
%	options.OUT =  1 print the observed p-values, 0 don't show it. OUT = 1 by
%                  default.
%   options.Pobs = 1 compute p-values only for observed (not permuted)
%                  data. This saves a considerable amount of time if the number
%                  of variables or permutations is high. Pobs = 0 by default
%                  (compute p-values for both observed and permuted data).
%   options.w =    1-by-m vector of weights for the combining function.
%                  W = [1 1 ...].*1/m by default. 
%   options.threshold = is a number [0;1]. The combining function does not 
%                  consider threshold values greater than OPTIONS.threshold
%                  (smaller if we are using 'D' or 'M' combining function).
%   options.labels.dims{i} = {'Var1', ......,'Varm'} to customize the
%                  labels of p-values in the output for the variables in
%                  dimension i of the T matrix. By default the labels are
%                  Y1, ..., Ym. With i = ndims(T)+1, we specify the labels
%                  of p-values in the output for the combined variable.
%       options.Combdims dimension of the combination. By default it is set
%                equal to ndims(T), i.e. the last dimension of 
%                matrix T.
%
%[P_sec,T_sec] = NPC(...) returns the matrix P_sec and the value of the
%combined test statistics in the (B+1) length matrix T_sec with one
%dimension less than T. The last row of T_sec is related to the observed
%sample.
%
%[...,OPTIONS] = NPC(...,OPTIONS) saves the options used for the analysis.
%OPTIONS is a structural array with the following structure:
%
%       options.labels.dimslabel label for the dimensions of the P-value
%                 matrix P:
%                 First dimension label 'Random Permutation'
%                 Second dimension label 'VARIABLES Y'
%                 Third dimension label 'VARIABLES X'
%                 Fourth dimension label 'STRATA'
%       options.labels.dims labels for any variable of any dimension of the
%                P-value matrix P
%       options.OUT see above
%       options.Pobs see above
%       options.Combdims dimension of the combination. By default it is set
%                equal to ndims(P), i.e. the last dimension of 
%                matrix P.
%       options.p.raw    P-values
%
%Example:
%We consider the data set GERMINA in the matlab file germinaNPC.mat.
%Below is the data set with two samples (FERTILIZER = 1 and FERTILIZER = 2)
%and four variables (Germinated, Weight, Surface and Surface2):
%
%           Fertilizer	Germinated	Weight	Surface	Surface2
%               0           1       3,88	8,09	65,4481
%               0           1       2,04	5,76	33,1776
%               0           1       5,48	18,01	324,3601
%               0           1       2,31	8,81	77,6161
%               0           1       1,9     8,17	66,7489
%               0           1       1,75	6,62	43,8244
%               0           1       3,02	7,69	59,1361
%               0           0			
%               0           0			
%               0           0			
%               1           1       3,93	19,29	372,1041
%               1           1       2,56	10,77	115,9929
%               1           1       8,3     18,81	353,8161
%               1           1       4,21	10,56	111,5136
%               1           1       1,86	9,48	89,8704
%               1           1       3,09	12,54	157,2516
%               1           1       5,09	18,35	336,7225
%               1           1       4,08	11,84	140,1856
%               1           1       3,63	11,44	130,8736
%               1           1       2,61	7,66	58,6756
%
%First we perform a two-samples test with the NP_2s function:
%
%load germinaNPC;
%y = germinaNPC(:,2:5);
%x = germinaNPC(:,1);
%options.labels.dims{3} = {'Fertilizer'};
%options.labels.dims{2} = {'Germinate','Weight','Surface','Surface2'};
%[P T] = NP_2s(y,x,1000,0,options);          
%
%which produce the output:
%  
%_______________________________________________________________
% Testing equality in distribution for 2 independent samples 
% for continuous (or dichotomous) variables (with or without missing values) 
% Number of Conditional Montecarlo: 1001
%            			            
%	VARIABLES Y 	Germinate   	Weight      	Surface     	Surface2    
%	  Fertilizer	0.2208      	0.2368      	0.07692     	0.1269      
%
%Now we can combine the four partial tests with NPC to obtain a global test
%using Fisher’s combining function:

%
%[P_sec, T_sec] = NPC(P,'F',options);
%
%which produce the output:
%
%________________________________________________________________
%Non Parametric Combination procedure
%Number of Conditional Montecarlo: 1001
%            			            
%	VARIABLES Y 	Fertilizer   
%	Germinate   	0.21578      
%	Weight      	0.21379      
%	Surface     	0.067932     
%	Surface2    	0.12188      
%
% Comb Funct		 Fisher
% p-GLOBAL	
%	            	Fertilizer  
%	            	0.074925    
%
%This function is part of NPClib.
%NPC Lib Matlab library (companion software of F. Pesarin, L. Salmaso: Permutation Tests for Complex Problems: theory, applications and software, John Wiley and sons).
%developed by Livio Finos; consulting team: Rosa Arboretti, Francesco Bertoluzzo, Stefano Bonnini, Chiara Brombin, Livio Corain, Luigi Salmaso and Aldo Solari.
%E-mail for technical questions: livio@stat.unipd.it (Livio Finos).

if nargin==1
    comb_funct='T';
end

if nargin==3
    options=get_options(T,'NPC',options);
else
    options=get_options(T,'NPC');
end

%%%%%%%%%%%%%%%nuovo
if iscell(comb_funct)
    options.p.raw=t2p(T,1);
else
    switch comb_funct
        case {'Direct','D','Max-t','M'}
            options.p.raw=t2p(T,1);
        otherwise
            options.p.raw=T(end,:);
    end
end
%%%%%%%%%%%%%%%%%%%
dims=size(T);
temp=ones(1,options.Combdims);
temp(1:length(dims))=dims;
dims=temp;
if isfield(options,'Combdims')
    T=split2dims2(T,options.Combdims);
    T=reshape(T,[dims(1) dims(options.Combdims) prod(dims(setdiff(2:length(dims),options.Combdims)))]);
else
    T=reshape(T,[dims(1:2) prod(dims(3:length(dims)))]);
end

B=size(T,1);

if iscell(comb_funct)

t=ones(1,0);
for i=1:length(comb_funct)
    if size(comb_funct{i},2)==1,
        CFF=cell(size(comb_funct{i}));
        CFF(:,1)=comb_funct{i}(:,1);
        comb_funct{i}=CFF;
    end
    for j=size(comb_funct{i},1):-1:1
        switch comb_funct{i}{j,1}
            case {'D'}
                comb_funct{i}{j,1}='Direct';
            case {'M'}
                comb_funct{i}{j,1}='Max-t';
            case {'F'}
                comb_funct{i}{j,1}='Fisher';
            case {'L'}
                comb_funct{i}{j,1}='Liptak';
            case {'LL'}
                comb_funct{i}{j,1}='LiptakLog';
            case {'Ma'}
                comb_funct{i}{j,1}='Mahalanobis';
            case {'T'}
                comb_funct{i}{j,1}='Tippett';
        end
        o.OUT=0;
        o.Pobs=0;
        if not(isempty(comb_funct{i}{j,2}))
            o.w=comb_funct{i}{j,2}(1:size(T,2));
            comb_funct{i}{j,2}=comb_funct{i}{j,2}(1:size(T,2));
        else 
            comb_funct{i}{j,2}=ones(1,size(T,2));
        end
        [T_sec(:,j,:) t_sec(:,j,:)]=NPC(T,comb_funct{i}{j,1},o);
    end
    t=[t;reshape(T_sec(end,:,:),prod(size(T_sec(end,:,:))),1)];
    T=T_sec;
    T2=t_sec;
    clear T_sec t_sec
end
P_sec=T;
T_sec=T2;
CF='Iterated';

else
    
switch comb_funct
    case {'Direct','D','Max-t','M'}
        if isfield(options,'threshold'),T(T<options.threshold)=0;end
    otherwise
        if isfield(options,'threshold'),T(T>options.threshold)=1;end
end

if isfield(options,'w'), if size(T,3)==size(options.w,3); end, end

switch comb_funct
    case {'Direct','D','d'}
        CF='Direct';
        if isfield(options,'w'),T_sec=sum(T.*repmat(options.w,[B 1 1]),2);
        else         T_sec=sum(T,2);end
       % options.p.raw=t2p(T,1);
       % options.p.raw=options.p.raw(:)';
    case {'Max-t','M','m'}
        CF='Max-t';
        if isfield(options,'w'),T_sec=max(T.*repmat(options.w,[B 1 1]),[],2);
        else T_sec=max(T,[],2);     end
        %options.p.raw=t2p(T,1);
        %options.p.raw=options.p.raw(:)';
    case {'Fisher','F','f'}
        CF='Fisher';
        %options.p.raw=T(end,:,:);
        T=(T.*B -.5)/B;
        if isfield(options,'w'),T_sec=-sum(log(T).*repmat(options.w,[B 1 1]),2);
        else T_sec=-sum(log(T),2);end
    case {'Liptak','L','l'}
        CF='Liptak';
        %options.p.raw=T(end,:,:);
        T=(T.*B -.5)/B;
        if isfield(options,'w'),T_sec=-sum(norminv(T).*repmat(options.w,[B 1 1]),2);
        else T_sec=-sum(norminv(T),2);end
    case {'LiptakLog','LL','ll'}
        CF='LiptakLog';
        %options.p.raw=T(end,:,:);
        T=(T.*B -.5)/B;
        if isfield(options,'w'),T_sec=sum(log((1-T)./T).*repmat(options.w,[B 1 1]),2);
        else T_sec=sum(log((1-T)./T),2);end
    case {'Mahalanobis','Ma','ma'}
        CF='Mahalanobis';
        %options.p.raw=T(end,:,:);
        if isfield(options,'w'),fprintf('\n***WARNING: weights are not consider for Mahalanobis Combining Function');end
        T=(T.*B -.5)/B;
        T=norminv(T);
        for idim=1:size(T,3)
        invcov=(T(:,:,idim)'*T(:,:,idim)./B)^-1;
        for i=B:-1:1
            T_sec(i,1,idim)=T(i,:,idim)*invcov*T(i,:,idim)';
        end
        end
        %T_sec=T_sec';
    case {'Tippett','T','t'}
        CF='Tippett';
        %options.p.raw=T(end,:,:);
        if isfield(options,'w'),T_sec=max(-T./repmat(options.w,[B 1 1]),[],2);
        else T_sec=-min(T,[],2); end
end

[P_sec T_sec]=t2p(T_sec,options.Pobs);
end

options.p.raw=reshape(options.p.raw,[1 dims(2:end)]);
if isfield(options,'Combdims')
    T_sec=reshape(T_sec,[dims(1:options.Combdims-1) 1 dims(options.Combdims+1:end)]);
    if options.Pobs==1
        dims(1)=1;
    end
    P_sec=reshape(P_sec,[dims(1:options.Combdims-1) 1 dims(options.Combdims+1:end)]);
    options.p.glb=reshape(P_sec(end,:,:),[1 dims(2:options.Combdims-1) 1 dims(options.Combdims+1:end)]);
%    options.labels.dims{options.Combdims}=options.labels.dimslabel{options.Combdims};
else
    T_sec=reshape(T_sec,[dims(1) 1 dims(3:end)]);
    options.p.glb=reshape(P_sec(end,:),[1 1 dims(3:end)]);
%    options.labels.dims{2}=options.labels.dimslabel(2);
end
    options.p.adj=[];
    options.B=B;


if options.OUT(1)==1        
    options.title='Non Parametric Combination procedure';
    options.CF=CF;
%    options.labels.dims{options.Combdims}=options.labels.dimslabel{options.Combdims};
    NPC_out(options);
    if iscell(comb_funct),
        g=1;
        fprintf('\nIterated combining procedure:\nIteration\tCombination\tCombining funct.\tP-value\t\tWeights')
    for i=1:length(comb_funct)
        for j=1:size(comb_funct{i},1)
                 fprintf('\n\t%2.0f\t\t%2.0f',i,j)
            fprintf('\t\t %10s\t\t\t\t%1.4f\t\t',comb_funct{i}{j,1},t(g))
            fprintf('%1.3f ',comb_funct{i}{j,2})
            g=g+1;
        end
        fprintf('\n')
    end     
    end
end

if nargout==3,  options.labels.dims{options.Combdims}=options.labels.dimslabel{options.Combdims}; end
