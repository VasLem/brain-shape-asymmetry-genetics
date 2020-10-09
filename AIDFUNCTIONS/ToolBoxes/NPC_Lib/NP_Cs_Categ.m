function [P, T, options, P_sub, T_sub] = NP_Cs_Categ(Y,X,B,var_type,comb_funct,options)
%Nonparametric permutation test for C > 2 samples for categorical (ordinal
%or not) variables.
%P = NP_Cs_Categ(Y,X) performs a one-way ANOVA to compare the
%distributions of C > 2 samples. Y is the n-by-m data matrix 
%(n=n1+n2+ ... +nC sample size, m variables number) of categorical
%responses. X is a vector of length n. The C modality of X indicates to which
%sample the corresponding row of Y belongs. If X is an n-by-q matrix, the same
%analysis is replicated for each column. If X or Y are cell arrays, the cells
%are labels of the variables of X and Y and variables must be declared global
%variables using reminD.m (and have a structure as in xlsimport.m, see example
%below).
%P is the (B+1)-by-m-by-q P-values matrix, where B is the number of
%random permutations set equal to 1000 by default. The last row of P is
%related to the observed sample. The alternative of the test is
%"distributions are not equal". 
%
%P = NP_Cs_Categ(Y,X,B) performs B random permutations. If B is a scalar,
%it indicates the number of random permutations. If B is a (B+1)-by-n matrix,
%it indicates the permutation sample space; in this case, each row is a
%(random) permutation of the vector 1:n. The last row has to be the vector
%1:n itself if we wish to make inference on the observed data. B = 1000
%by default.  
%
%P = NP_Cs_Categ(Y,X,B,VAR_TYPE) to specify the kind of categorical
%variables. VAR_TYPE is set equal to 1 for non ordinal variables, and equal to
%2 for ordinal variables. If VAR_TYPE is a scalar, the same type is assumed
%for each variable. If it is a vector of length m, different types are 
%specified for each variable. By default VAR_TYPE is set at 1.
%
%P = NP_Cs_Categ(Y,X,B,VAR_TYPE,COMB_FUNCT) to specify the combining
%function to use to combine the partial tests. Possible settings for
%COMB_FUNCT are:
%
%           'D' for Direct combination of t-statistic,
%           'M' for Max-t  combining function of t-statistic,
%           'F' for Fisher's combining function of p-values,
%           'L' for Liptak's combining function of p-values,
%           'T' for Tippett's (min-p) combining function of p-values.
%           Extras:
%           'AD' performs the Anderson-Darling test for ordinal variables. 
%                For tail==0 'AD' is based on the absolute values of the
%                statistics. If we wish to base on the square values,
%                write 'A2'.
%           'KS' performs the Kolmogorov-Smirnov test.
%           'X2' performs the Chi square (exact) test.
%
%By default COMB_FUNCT is set at 'F'.
%
%P = NP_Cs_Categ(...,'OPTIONS') Possible settings for options are:
%
%   options.OUT is a 1-by-2 vector. The first element is related to the higher
%                  level analysis, the second one is related to the lower level
%                  analysis (p-values related to the categories). If 1, print the
%                  observed p-values, if 0 don't show it. OUT = [1 0] by
%                  default.
%   options.Pobs = 1 compute p-values only for observed (not permuted)
%                  data. This saves a considerable amount of time if the number
%                  of variables or permutations is high. Pobs = 0 by default
%                  (compute p-values for both observed and permuted data).
%   options.labels.dims{2} = {'Var1', ......,'Varm'} to customize the
%                  labels of p-values in the output for the m variables Y.
%                  By default the labels are Y1, ..., Ym
%   options.labels.dims{3} = {'Var1', ......,'Varq'} to customize the
%                  labels of p-values in the output for the q variables X.
%                  By default the labels are X1, ..., Xq
%   options.tail   specify the alternative hypotheses when there are 2 
%                  samples. The null hypothesis is: "distributions are equal". 
%                  Possible settings for options.tail are:
%
%                  0,  "distributions are not equal" (always when variables
%                       are nominal)
%                  1,  "CDF of first group is less than CDF of the second group"
%                 -1,  "CDF of first group is greater than CDF of the second group"
%   
%                  By default options.tail is set at 0. If options.tail is a
%                  vector (or a matrix if X is a matrix), each element
%                  corresponds
%                  to a variable; if is a scalar, this is considered for each
%                  variable.
%
%[P,T] = NP_Cs_Categ(...) returns the P matrix and value of the test statistics
%in the (B+1)-by-m-by-q T matrix. The last row is related to the observed
%sample.
%
%[...,OPTIONS] =  NP_Cs_Categ(...) saves the options used for the analysis. 
%OPTIONS is a structural array with the following structure:
%
%       options.labels.dimslabel label for the dimensions of the P-value
%                matrix P:
%		         First dimension label 'Random Permutation'
%		         Second dimension label 'VARIABLES Y'
%		         Third dimension label 'VARIABLES X'
%		         Fourth dimension label 'STRATA'
%	    options.labels.dims labels for any variable of any dimension of the
%                P-value matrix P
%       options.OUT see above
%       options.Pobs see above
%       options.Combdims dimension of the combination. By default it is set
%                equal to length(size(P)), i.e. the last dimension of 
%                matrix P
%	    options.p.raw    P-values
%
%[...,P_sub] =  NP_Cs_Categ(...) also returns the (B+1)x(C1+C2+...+Cm) 
%p-values matrix, with C total number of dummy variables from the original
%number m for the multi-focus (independence) analysis (see references).
%The last row is related to the observed sample.
%
%[...,P_sub,T_sub] =  NP_Cs_Categ(...) to have the statistics matrix of the 
%multi-focus analysis. The last row is related to the observed one.
%
%The permutation methodology can handle missing data 
%(whether missing at random or not). Please indicate this as NaN.
%
%Example:
%We consider the data set PERCRAWDATA in the matlab file
%percrawdataCsCateg.mat. Below is the data set with four samples (identified by
%variable dose which has values 0, 150, 1500, 5000) and four response
%categorical variables (hc_posture, removal, handling, lacrimation). The
%total number of subject is 25:
%
%   rat     dose	hc_posture	removal	handling	lacrimation
%   1       150     6           1       2           2
%   2       1500	2           2       4           2
%   3       150     3           1       2           1
%   4       150     3           2       2           1
%   5       0       1           1       2           1
%   6       5000	4           1       1           2
%   7       0       3           1       2           1
%   8       150     3           1       2           1
%   9       1500	2           2       2           1
%   10      500     1           1       2           1
%   11      1500	2           2       3           1
%   12      0       1           1       2           1
%   13      5000	5           1       1           3
%   14      150     3           2       3           1
%   15      1500	3           1       2           1
%   16      150     3           1       2           1
%   17      1500	3           1       2           2
%   18      0       1           1       2           1
%   19      0       3           1       2           1
%   20      0       1           1       2           1
%   21      500     1           2       2           1
%   22      150     2           2       2           2
%   23      500     2           2       4           2
%   24      500     1           1       2           1
%   25      1500	1           1       2           1
%
%To perform a C-sample test with these data, the following vectors are
%used as input:
%
%load percrawdataCsCateg;
%x = percrawdataCsCateg(:,2);
%y = percrawdataCsCateg(:,3:end);
%[P,T,options] = NP_Cs_Categ(y,x);
%
%which produce the output:
%
% Testing C independent samples for categorical variables: p-values 
%			 VARIABLES Y
%           		Y1        	Y2        	Y3        	Y4        
%           		0.007992   	0.24575    	0.037962   	0.098901
%
%To customize the labels of X and Y variables, we must define the
%variables labels in the options:
%
%options.labels.dims{2} = {'hc_posture','removal','handling','lacrimation'};
%options.labels.dims{3} = {'dose'};
%
%and repeat the test with the options (if options is used, it must also
%specify the B,var_type and comb_funct parameters);
%
%[P,T] = NP_Cs_categ(y,x,B,var_type,comb_funct,options);
%
%which produce the output:
%
%Testing C independent samples for categorical variables: p-values 
%			 VARIABLES Y
%           		hc_posture	removal   	handling  	lacrimation        
%      dose			0.007992   	0.24575    	0.037962    0.098901
%
%These data are also reported in the Excel file percrawdataCsCateg.xls,
%with column labels showing the names of the variables. To perform the same
%analysis with this file we use the structural arrays D generated with
%xlsimport:
%
%[D] = xlsimport('percrawdataCsCateg');
%[P,T,options] = NP_Cs_Categ(D(:,3:end),D(:,2));
%
%This produces the output:
%
%               VARIABLES Y
%                   hc_posture	removal   	handling  	lacrimation
%      dose         0.007992   	0.24575    	0.037962    0.098901  
%
%It is possible to perform the same analysis, whith the same output as above,
%using the variable labels of the columns of the Excel files as input.
%Beforehand, to obtain this result, it is necessary to declare the variable
%D global:
%
%[D]=xlsimport('percrawdataCsCateg');
%reminD(D)
%[P,T,options] = NP_cs({'hc_posture','removal','handling', ...
%   'lacrimation'},'dose');
%
%This function is part of NPClib.
%NPC Lib Matlab library (companion software of F. Pesarin, L. Salmaso: Permutation Tests for Complex Problems: theory, applications and software, John Wiley and sons).
%developed by Livio Finos; consulting team: Rosa Arboretti, Francesco Bertoluzzo, Stefano Bonnini, Chiara Brombin, Livio Corain, Luigi Salmaso and Aldo Solari.
%E-mail for technical questions: livio@stat.unipd.it (Livio Finos).

if nargin == 2
   B=1000;
   var_type=1;
   comb_funct='F';
elseif nargin == 3
   var_type=1;
   comb_funct='F';
elseif nargin == 4
   comb_funct='F';
end


if nargin==6
    [Y, options]=getDopts(Y,2,options);
else
    [Y, options]=getDopts(Y,2);
end
[X, options]=getDopts(X,3,options);
sizY=size(Y,2);
sizX=size(X,2);
options=get_options(ones([1 sizY sizX]),'partial',options);

if isscalar(options.tail)
    options.tail=repmat(options.tail,[1 size(Y,2) size(X,2)]);
elseif isvector(options.tail)
    options.tail=repmat(options.tail(:)',[1 1 prod(sizX(2:end))]);
end


if length(var_type)==1
    var_type=var_type.*ones(1,size(Y,2));
end

[N m]=size(Y);
if  isscalar(B)
    [B] = Space_perm(N,B);
end


if (sum(abs(comb_funct-'X2')))
    [Y, orig,w,labels] = expand_categ(Y,var_type,options.labels.dims{2});
    options.orig=orig;
    if not(isfield(options,'tail_categ'))|(isfield(options,'tail'))
        n_categ=unique(orig)';
        options.tail_categ=ones(1,0,size(options.tail,3));
        for i=1:length(n_categ)
           options.tail_categ=[options.tail_categ   repmat(options.tail(1,i,:),[1 sum(orig-n_categ(i)==0) 1])];
        end
    end
end


if comb_funct=='X2'
    if    sum(isnan(Y(:)))==0
        [P T_sub no P_sub]  = NP_Chi2(Y,X,B,options);
    else
        fprintf('\n WARNING: missing values are not allowed with exact X2 test, choose any other comb_funct\n')
        P=NaN.*ones(1,size(Y,2));
        T_sub=P;
        P_sub=P;
    end
else
    opts=options;
    opts.Combdims=2;
    if isfield(opts,'tail_categ') opts.tail=opts.tail_categ; end
    opts.OUT=0;
    if (comb_funct(1)=='A'| prod((comb_funct=='KS')*1)| comb_funct=='D' | comb_funct=='M'), opts.Pobs=1; end
    if length(unique(X))==2    
        [P_sub T_sub]=NP_2s(Y,X,B,options.tail_categ,opts);
    else
        [P_sub T_sub]=NP_Cs(Y,X,B,opts);
    end 
end
opts.OUT=options.OUT;

if length(opts.OUT)==1
    opts.OUT=0;
else
    opts.OUT=options.OUT(2);
end
opts.Pobs=options.Pobs;
P=ones(size(B,1), length(options.labels.dims{2}),size(X,2));
T=zeros(size(B,1), length(options.labels.dims{2}),size(X,2));
if comb_funct(1)=='A'
    for i=unique(orig)
        opts.labels.dims{2}=labels(orig==i);
        opts.labels.dimslabel{2}=options.labels.dims{2}{1};
        if opts.OUT==1, fprintf(['\n --- Partial p-value for variable: ' options.labels.dims{2}{i} '\n']);end
        N=nansum(Y(:,orig==i));
        opts.w=[N.*(size(Y,1)-sum(isnan(Y(:,orig==i))) -N)].^.5;
        if prod((comb_funct=='AD')*1)
            [P(:,i,:) T(:,i,:)]=NPC(T_sub(:,orig==i).^2,'D',opts);
        else
            [P(:,i,:) T(:,i,:)]=NPC(T_sub(:,orig==i),'D',opts);
        end
    end
elseif comb_funct=='KS'
    for i=unique(orig)
        opts.w=ones(size(T_sub(1,orig==i)));
        opts.labels.dims{2}=labels(orig==i);
        opts.labels.dimslabel{2}=options.labels.dims{2}{1};
        if opts.OUT==1, fprintf(['\n --- Partial p-value for variable: ' options.labels.dims{2}{i} '\n']);end
        
        [P(:,i,:) T(:,i,:)]=NPC(abs(T_sub(:,orig==i)),'M',opts);
    end
elseif comb_funct=='X2',
    
elseif comb_funct=='D' | comb_funct=='M'
    for i=unique(orig)
        opts.labels.dims{2}=labels(orig==i);
        opts.labels.dimslabel{2}=options.labels.dims{2}{1};
        opts.w=ones(size(T_sub(1,orig==i)));
        if opts.OUT==1, fprintf(['\n --- Partial p-value for variable: ' options.labels.dims{2}{i} '\n']);end
        [P(:,i,:) T(:,i,:)]=NPC(T_sub(:,orig==i,:),comb_funct,opts);
    end
else
    for i=unique(orig)
        opts.w=ones(size(P_sub(1,orig==i,:)));
        opts.labels.dims{2}=labels(orig==i);
        opts.labels.dimslabel{2}=options.labels.dims{2}{1};
        if opts.OUT==1, fprintf(['\n --- Partial p-value for variable: ' options.labels.dims{2}{i} '\n']);end
        [P(:,i,:) T(:,i,:)]=NPC(P_sub(:,orig==i,:),comb_funct,opts);
    end
end

P(:,setdiff(1:size(P,2),unique(orig)),:)=1;
T(:,setdiff(1:size(P,2),unique(orig)),:)=0;

options.p.raw=P(end,:,:);
NP_out('Testing C indipendent samples for categorical variables',options)
