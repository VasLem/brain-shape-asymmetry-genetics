function [P,T, options,P_sub,T_sub,pp,tt] = NP_StOrd(Y,X,B,tail,comb_funct,options)
%Nonparametric permutation test for stochastic ordering. 
%P = NP_StOrd(Y,X) performs m nonparametric permutation tests of the null
%hypothesis of absence of stochastic ordering against the alternative
%hypothesis of monotonic (increasing or decreasing) effect on the C > 1
%independent samples. Y is the n-by-m (n=n1+n2+...nC) data matrix (n sample
%size, m variables number) of responses. X is a vector of length n*C.
%The modality of X indicates to which sample the corresponding row of Y belongs.
%P is the (B+1)-by-m P-values matrix, where B is the number of random
%permutations set equal to 1000 by default. The last row of P is
%related to the observed sample. The default alternative of the test is
%"there is an increasing effect". We can change the default alternative by
%parameter TAIL (see below).
%
%P = NP_StOrd(Y,X,B) performs B random permutations. If B is a scalar,
%it indicates the number of random permutations. If B is a (B+1)-by-n matrix,
%it indicates the permutation sample space; in this case, each row is a
%(random) vector of 1 and -1. The last row must be the vector ones(1,n).
%If B is equal to 0, all possible permutations are done.
%
%P = NP_StOr(Y,X,B,TAIL) The default alternative is "increasing effect". 
%With TAIL it is possible to specify other alternatives: 
%   
%	     TAIL =  0, alternative: "means are not equal"
%        TAIL =  1, alternative: "increasing effect"
%        TAIL = -1, alternative: "decreasing effect"
%
%If TAIL is a vector (or a matrix if X is a matrix), each element
%corresponds to a variable; if it is a scalar, this is considered for each
%variable.
%
%P = NP_StOr(Y,X,B,comb_funct) to specify the combining
%function to use to combine the partial tests. Possible settings for
%COMB_FUNCT are:
%
%           'D' for Direct combination of t-statistic,
%           'M' for Max-t combining function of t-statistic,
%           'F' for Fisher's combining function of p-values,
%           'L' for Liptak's combining function of p-values,
%           'T' for Tippett's (min-p) combining function of p-values.
%           extras:
%           'AD' performs the Anderson-Darling test for ordinal variables
%                for tail==0, 'AD' is based on the absolute values of the
%                statistics. If we wish to base it on the square values,
%                write 'A2'.
%           'KS' performs the Kolmogorov-Smirnov test.
%           'X2' performs the Chi square (exact) test.
%
%By default COMB_FUNCT is set to 'F'.
%
%P = NP_StOr(...,'OPTIONS') Possible settings for options are:
%
%   options.OUT is a 1-by-2 vector. The first element is related to the higher
%                  level analysis, the second one is related to the lower level
%                  analysis (p-values related to the categories). 1 print the
%                  observed p-values, 0 don't show it. OUT = [1 0] by
%                  default.
%   options.Pobs = 1 compute p-values only for observed (not permuted)
%                  data. This saves a considerable amount of time if the number
%                  of variables or permutations is high. Pobs = 0 by default
%                  (compute p-values for both observed and permuted data).
%   options.labels.dims{2} = {'Var1', ......,'Varm'} to customize the
%                  labels of p-values in the output for the m variables Y.
%                  By default the labels are Y1, ..., Ym
%   options.labels.dims{3} = {'Label first comparison', ...,'Label K-th
%                  comparison'} to customize the labels of p-values in the
%                  output for the K comparisons.
%   options.match  matching n-by-1 vector for related samples. Not matching
%                  as the default.
%   options.CF_match = Combining Function to use for the combination of
%                  separated tests of each matching. OPTIONS.CF_match='D'
%                  as the default.
%
%[P,T] = NP_StOrd(...) also returns the values of the test statistics
%in the (B+1)-by-m T matrix. The last row is related to the observed
%sample.
%
%[...,OPTIONS] = NP_StOrd(...) saves the options used for the analysis.
%OPTIONS is a structural array with the following structure:
%
%       options.labels.dimslabel label for the dimensions of the P-value
%                matrix P:
%                First dimension label 'Random Permutation'
%                Second dimension label 'VARIABLES Y'
%                Third dimension label 'VARIABLES X'
%                Fourth dimension label 'STRATA'
%       options.labels.dims labels for any variable of any dimension of the
%                P-value matrix P
%       options.OUT see above
%       options.Pobs see above
%       options.Combdims dimension of the combination. By default it is set
%                equal to length(size(P)), i.e. the last dimension of 
%                matrix P
%       options.p.raw    P-values
%
%[...,P_sub] = NP_StOrd(...) also return the (B+1)-by-(C-1)-by-m p-values
%matrix, with C total number of dummy variables from the original number
%of m for the multi-focus (independence) analysis (see references). The
%last row is related to the observed sample.
%
%[...,P_sub,T_sub] = NP_StOrd(...) to have the statistics matrix of the
%multi-focus analysis. The last row is related to the observed sample.
%
%The permutation methodology can handle missing data (whether missing at
%random or not). Please indicate this as NaN.
%
%Example:
%We consider the data set PERCRAWDATA in the matlab file
%percrawdataStOrd.mat. Below is the data-set with four samples (identified by
%variable dose which assumes the values 0, 150, 1500, 5000) and four response
%variables ('foot_splay','forelimb_grip','hindlimb_grip','weight'). The total
%number of subjects is 24:
%
%       rat	dose	foot_splay	forelimb_grip	hindlimb_grip	weight
%       1	150     20,5        0,43            0,2             140,7
%       2	1500	35          0,81            0,245           144,6
%       3	150     25,5        0,86            0,295           151,7
%       4	150     22          0,65            0,365           161,1
%       5	0       31          0,71            0,295           144,8
%       6	0       32          0,425           0,415           154
%       7	150     34          0,715           0,395           152,1
%       8	1500	26          0,475           0,215           147
%       9	500     39,5        0,84            0,39            161,2
%       10	1500	39,5        0,525           0,195           147,6
%       11	0       21          0,49            0,34            153,4
%       12	150     33          0,635           0,355           151,5
%       13	5000	33          0,55            0,25            154,1
%       14	1500	35,5        0,575           0,255           149,9
%       15	150     39          0,84            0,285           160,4
%       16	1500	36          0,635           0,2             139,2
%       17	0       34          0,595           0,15            148,6
%       18	0       29          0,58            0,22            158,8
%       19	0       32,5        0,635           0,195           151
%       20	500     36,5        0,76            0,22            144,5
%       21	150     36,5        0,75            0,36            150,9
%       22	500     34          0,79            0,23            142,2
%       23	500     38,5        0,465           0,215           151
%       24	1500	27,5        0,745           0,11            152
%
%To test the presence of stochastic ordering:
%
%load percrawdataStatOrd
%y = percrawdataStatOrd(:,3:end);
%x = percrawdataStatOrd(:,2);
%
%which produce the output:
%
%_______________________________________________________________
%
%Testing Stochastic Ordering 
%			        VARIABLES Y
%           		Y1        	Y2        	Y3        	Y4        
%           		0.84316    	0.5005     	0.064935   	0.17582  
%
%to customize the labels of Y variables, we must define the variable
%labels in the options:
%
%options.labels.dims{2} = {'foot_splay','forelimb_grip','hindlimb_grip', ...
%                    'weight'};
%
%and repeat the test with the options (if options is used, it must
%also specify the B, TAIL and the comb_funct parameters):
%
%[P,T,options] = NP_StOrd(y,x,1000,1,'F',options); 
%
%which produce the output:
%________________________________________________________________
%
% Testing Stochastic Ordering 
%               VARIABLES Y
%          		foot_splay	fore.._grip	hind.._grip	weight    
%           		0.85215    	0.53846    	0.041958   	0.15385   
%
%These data are also reported in the Excel file percrawdataStOrd.xls, with
%column labels showing the names of the variables. To perform the same
%analysis whit this file, we use the structural arrays D generated with
%xlsimport:
%
%[D] = xlsimport('percrawdataStOrd');
%[P,T,options] = NP_StOrd(D(:,3:end),D(:,2));
%
%This produce the output:
%________________________________________________________________
%
% Testing Stochastic Ordering 
%                   VARIABLES Y
%           		foot_splay	fore.._grip	hind.._grip     weight    
%      dose         0.88312    	0.54146    	0.044955        0.16583
%
%It is possible to perform the same analysis, with the same output as above,
%using the variable labels of the columns in the Excel files as input.
%Beforehand, it is necessary to declare the variable D global:
%
%[D]=xlsimport('PERCRAWDATAStatOrd');
%reminD(D)
%[P,T,options] = NP_cs({'foot_splay','forelimb_grip','hindlimb_grip', ...
%            'weight'},'dose');
%
%This function is part of NPClib.
%NPC Lib Matlab library (companion software of F. Pesarin, L. Salmaso: Permutation Tests for Complex Problems: theory, applications and software, John Wiley and sons).
%developed by Livio Finos; consulting team: Rosa Arboretti, Francesco Bertoluzzo, Stefano Bonnini, Chiara Brombin, Livio Corain, Luigi Salmaso and Aldo Solari.
%E-mail for technical questions: livio@stat.unipd.it (Livio Finos).

if nargin==2
    B=1000;
    tail=0;
    comb_funct='D';
elseif nargin==3
    tail=0;
    comb_funct='D';
elseif nargin==4
    comb_funct='D';
end

if nargin==6
    [Y, options]=getDopts(Y,2,options);
else
    [Y, options]=getDopts(Y,2);
end
[X, options]=getDopts(X,3,options);

    sizY=size(Y,2);
    sizX=size(X,2);
    options=get_options(ones(1,sizY,sizX),'StOrd',options);



if nargin <= 1
    error('Not enough iNPut argument')
elseif nargin == 2
   B=1000;
   tail=1;
   comb_funct='F';
elseif nargin == 3
   tail=1;
   comb_funct='F';
elseif nargin == 4
   comb_funct='F';
end

if sum(tail==0)
    fprintf('\nWarning: one or more tails are set to 0.\nThis is not allowed in Stochastic ordering, these are now set to 1 (i.e. positive trend).\n')
end
if isscalar(B)
    perm = Space_perm(size(Y,1),B);
    B=B+1;
else
    perm=B;
    B=size(B,1);
end


if not(isfield(options,'matchs'))
    options.matchs=ones(size(Y,1),1);
else
    if (isempty(options.matchs))
        options.matchs=ones(size(Y,1),1);
    end
end



[X orig]=expand_categ(X,2,options.labels.dims{3});

matchs=unique(options.matchs);

P_sub=zeros(B,size(Y,2),size(X,2),length(matchs));
T_sub=zeros(B,size(Y,2),size(X,2),length(matchs));
for i=length(matchs):-1:1
    if length(unique(X(find(options.matchs==matchs(i)),:)))==2
        [P_sub(:,:,:,i), T_sub(:,:,:,i)] = NP_2s(Y(find(options.matchs==matchs(i)),:),X(find(options.matchs==matchs(i)),:),perm,tail,0);
    else
        P_sub(:,:,:,i)=1; T_sub(:,:,:,i)=0;
    end
end
pp=P_sub;
tt=T_sub;

if length(unique(options.matchs))>1
opts.OUT=0;
opts.Combdims=4;
switch options.CF_match
    case {'Direct','D','Max-t','M'}
        if strmatch(comb_funct,{'Direct','D','Max-t','M'}), opts.Pobs=1; end
        [P_sub T_sub]=NPC(T_sub,options.CF_match,opts);
    case {'Fisher','F','Liptak','L','Tippett','T'}
        [P_sub T_sub]=NPC(P_sub,options.CF_match,opts);
end
end


opts.OUT=0;
opts.Combdims=3;
opts.Pobs=options.Pobs;
for i=length(unique(orig)):-1:1
    switch comb_funct
        case {'Direct','D','Max-t','M'}
            [P(:,:,i) T(:,:,i)] =NPC(T_sub(:,:,find(orig==i)),comb_funct,opts);
        case {'Fisher','F','Liptak','L','Tippett','T'}
            [P(:,:,i) T(:,:,i)]=NPC(P_sub(:,:,find(orig==i)),comb_funct,opts);
    end
end
  
options.p.raw=P(end,:,:);
NP_out('Testing Stochastic Ordering',options) 