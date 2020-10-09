function [P, T, options, P_sub, T_sub] = NP_ReM_Categ(Y,X,DES,B,var_type,comb_funct,options)
%Nonparametric permutation test for repeated measures of categorical
%(ordinal or not) variables.
%P = NP_ReM_Categ(Y,X,DES) performs m nonparametric permutation tests of the
%null hypothesis of equality of the mean of the m response categorical
%variables observed on a number of occasions, usually according to time,
%so that successive responses are dependent. Y is the [n*C]-by-m data
%matrix (n sample size, C number of repeated measures, m variables number)
%of categorical responses. X is a vector of length n*C. The modality of X
%indicates to which time the corresponding row of Y belongs. If X is a scalar,
%it indicates the number of repeated measures C. DES is the C-by-K
%Design Matrix of comparisons (K is the number of comparisons). By default
%DES is the C-by-C-1 matrix:
%
%                            [ 1  0  0 ...
%                             -1  1  0 ...
%                              0 -1  1 ...]
%
%where 1 and -1 in any column indicate which measurements to compare. With
%the default setting each variable at time t is compared with the same
%variable at time t+1. If DES is a string, the following values are allowed:
%
%                  All: All comparisons
%                  Bal: Base line vs others (Base line is the first one)
%                  Seq: Sequential comparisons (1vs2, 2vs3, ...C-1vsC)
%                  Trd: Trend (1vs2...C, 12vs3...C, 123vs4...C)
%
%by default DES = 'Seq'.
%P is the (B+1)-by-m-by-K P-values matrix, where B is the number of random
%permutations set at 1000 by default. The last row of P is
%related to the observed sample. The default alternative of the test is
%"means are not equal". We can change the default alternative by
%parameter options.TAIL (see below).
%
%P = NP_ReM_Categ(Y,X,DES,B) performs B random permutations.If B is a scalar,
%it indicates the number of random permutations. If B is a (B+1)-by-n matrix,
%it indicates the permutation sample space; in this case, each row is a
%(random) vector of 1 and -1. The last row must be the vector ones(1,n).
%If B is equal to 0 all possible permutations are done.
%
%P = NP_ReM_Categ(Y,X,DES,B,var_type) to specify the kind of categorical
%variables. VAR_TYPE is set equal to 1 for non ordinal variables, and equal to
%2 for ordinal variables. If it is a scalar, the same type is assumed
%for each variable, whereas if is a vector of length m, different types are
%specified for each variable. By default VAR_TYPE is set at 1.
%
%P = NP_ReM_Categ(Y,X,DES,B,var_type,comb_funct)to specify the combining
%function to use to combine the partial tests. Possible settings for
%COMB_FUNCT are:
%
%           'D' for Direct combination of t-statistic,
%           'M' for Max-t  combining function of t-statistic,
%           'F' for Fisher's combining function of p-values,
%           'L' for Liptak's combining function of p-values,
%           'T' for Tippett's (min-p) combining function of p-values.
%           extras:
%           'AD' performs the Anderson-Darling test for ordinal variables
%                for tail==0 'AD' is based on the absolute values of the
%                statistics.If we wish to base it on the square values,
%                write 'A2'. If missing data are present, these data are
%                not considered in the weight computations.
%           'KS' performs the Kolmogorov-Smirnov test. If missing data are
%                present, these data are not considered in the weight
%                computations.
%           'X2' performs the Chi square (exact) test.
%
%By default COMB_FUNCT is set to 'F'.
%
%P = NP_ReM_Categ(...,'OPTIONS') Possible settings for options are:
%
%   options.OUT is a 1-by-2 vector. The first element is related to the higher
%                  level analysis, the second to the lower level
%                  analysis (p-values related to the categories). 1 print the
%                  observed p-values, 0 don't show it. OUT = [1 0] by
%                  default.
%   options.Pobs = 1 compute p-values only for observed (not permuted)
%                  data. This save a considerable amount of time if the number
%                  of variables or permutations is high. Pobs = 0 by default
%                  (compute p-values for both observed and permuted data).
%   options.labels.dims{2} = {'Var1', ......,'Varm'} to customize the
%                  labels of p-values in the output for the m variables Y.
%                  By default the labels are Y1, ..., Ym
%   options.labels.dims{3} = {'Label first comparison', ...,'Label K-th
%                  comparison'} to customize the labels of p-values in the
%                  output for the K comparisons.
%   options.tail   specify the alternative hypotheses when there are 2
%                  samples. The null hypothesis is: "distributions are equal".
%                  Possible settings for options.tail are:
%
%                  0,  "distributions are not equal" (always when variables
%                       are nominal)
%                  1,  "CDF of first group is less than CDF of the second group"
%                 -1,  "CDF of first group is greater than CDF of the
%                       second group"
%
%                  By default options.tail is set at 0. If options.tail is a
%                  vector (or a matrix if X is a matrix), each element
%                  corresponds to a variable; if it is a scalar, this is considered
%                  for each variable.
%
%[P,T] = NP_ReM_Categ(...) also returns the values of the test statistics
%in the (B+1)-by-m-by-K T matrix. The last row is related to the observed
%sample.
%
%[...,OPTIONS] =  NP_Rem_Categ(...) saves the options used for the analysis.
%OPTIONS is a structural array with the following structure:
%
%       options.labels.dimslabel label for the dimensions of the P-value
%                matrix P:
%                 First dimension label 'Random Permutation'
%                 Second dimension label 'VARIABLES Y'
%                 Third dimension label 'VARIABLES X'
%                 Fourth dimension label 'STRATA'
%       options.labels.dims labels for any variable of any dimension of the
%                P-value matrix P
%       options.OUT see above
%       options.Pobs see above
%       options.Combdims dimension of the combination. By default it is set
%                equal to length(size(P)), i.e. the last dimension of 
%                matrix P
%       options.p.raw    P-values
%
%[...,P_sub] =  NP_Rem_Categ(...) also returns the (B+1)x(C1+C2+...+Cm)
%p-values matrix, with C total number of dummy variables from the original
%number of m for the multi-focus (independence) analysis (see references).
%The last row is related to the observed sample.
%
%[...,P_sub,T_sub] =  NP_Rem_Categ(...) to have the statistics matrix of the
%multi-focus analysis. The last row is related to the observed one.
%
%permutation methodology can handle missing data
%(whether missing at random or not). Please indicate this as NaN.
%
%Example:
%We consider the data set GUARDA in the matlab file guardaRemCateg.mat.
%Below the data set with 4 repeated measures in 7 subjects with 4 variables:
%
%               paz    Time    DM    DF    DR    LF
%               2    0       8    0    5    2
%               5    0       5    0    7    2
%               6    0       5    3    0    1
%               7    0       5    2    6    2
%               10    0       7    4    5    3
%               11    0       2    2    5    1
%               15    0       0    0    0    0
%               2    1       6    0    5    1
%               5    1       4    0    5    2
%               6    1       5    2    0    1
%               7    1       6    2    5    2
%               10    1       7    4    5    3
%               11    1       2    2    3    1
%               15    1       0    0    0    1
%               2    2       8    0    5    2
%               5    2       5    0    7    2
%               6    2       4    2    0    1
%               7    2       4    0    6    2
%               10    2       5    4    5    3
%               11    2       2    2    5    1
%               15    2       0    0    0    1
%               2    3       8    0    5    2
%               5    3       6    0    7    2
%               6    3       5    2    0    1
%               7    3       6    2    7    3
%               10    3       7    4    7    3
%               11    3       4    2    5    1
%               15    3       0    0    0    0
%
%To test the equality of means (means equal to 0) of all possible pairs of
%measurements and for each variable:
%
%load guardaRemCateg;
%x = guardaRemCateg(:,2);
%y = guardaRemCateg(:,3:end);
%[P,T,options] = NP_Rem_Categ(y,x);
%
%which produce the output:
%________________________________________________________________
%
%       Testing repeated measures for categorical variables: p-values
%                   VARIABLES
%       COMPARISONS        Y1            Y2            Y3            Y4
%           1vs2        0.25175        1              1              0.49151
%           2vs3        0.68232        1              1              1
%           3vs4        0.69031        1              0.4985         1
%           1vs3        0.49151        1              1              1
%           2vs4        1              1              0.25974        0.49151
%           1vs4        0.5045         1              0.4985         1
%
%to customize the labels of Y variables and the comparisons, we must
%define the variable %labels in the options:
%
%options.labels.dims{2} ={'DM','DF','DR','LF'};
%options.labels.dims{3} = {'T1 vs T2','T2 vs T3','T3 vs T4', ...
%    'T1 vs T3','T2 vs T4','T1 vs T4'};
%
%and repeat the test with the options (if options is used, it must
%also specify the B, TAIL and the comb_funct parameters):
%
%[P,T,options] = NP_Rem_categ(y,x,'All',1000,2,'F',options);
%
%which produce the output:
%
%Testing repeated measures for categorical variables: p-values
%________________________________________________________________
%
%       Testing repeated measures for categorical variables: p-values
%                    VARIABLES
%       COMPARISONS        DM            DF            DR            LF
%       T1 vs T2        0.25175        1              1              0.49151
%       T2 vs T3        0.68232        1              1              1
%       T3 vs T4        0.69031        1              0.4985         1
%       T1 vs T3        0.49151        1              1              1
%       T2 vs T4        1              1              0.25974        0.49151
%       T1 vs T4        0.5045         1              0.4985         1
%
%These data are also reported in the Excel file guardaRemCateg.xls, with column
%labels showing the names of the variables. To perform the same analysis 
%with this file, we use the structural arrays D generated with xlsimport:
%
%D = xlsimport('guardaRemCateg');
%reminD(D);
%options.labels.dims{3} = {'T1 vs T2','T2 vs T3','T3 vs T4', ...
%    'T1 vs T3','T2 vs T4','T1 vs T4'};
%[P,T,options] = NP_ReM({'DM','DF','DR','LF'},'Time','All');
%
%This function is part of NPClib.
%NPC Lib Matlab library (companion software of F. Pesarin, L. Salmaso: Permutation Tests for Complex Problems: theory, applications and software, John Wiley and sons).
%developed by Livio Finos; consulting team: Rosa Arboretti, Francesco Bertoluzzo, Stefano Bonnini, Chiara Brombin, Livio Corain, Luigi Salmaso and Aldo Solari.
%E-mail for technical questions: livio@stat.unipd.it (Livio Finos).

if nargin == 2
   B=1000;
   var_type=1;
   comb_funct='F';
   DES='Seq';
elseif nargin == 3
   B=1000;
   var_type=1;
   comb_funct='F';
elseif nargin == 4
   var_type=1;
   comb_funct='F';
end

if nargin==7
    [Y, options]=getDopts(Y,2,options);
else
    [Y, options]=getDopts(Y,2);
end
[X]=getDopts(X,3);
sizY=size(Y,2);
sizX=size(X,2);

if isempty(DES),      options.DES=DES_ReM_std(length(unique(X)),'All'); 
elseif ischar(DES),       options.DES=DES_ReM_std(length(unique(X)),DES);
else options.DES=DES;
end

options=get_options(ones(size(options.DES,1),size(Y,2),size(X,3)),'ReM',options);


if isscalar(options.tail)
    options.tail=repmat(options.tail,[1 size(Y,2) size(options.DES,2)]);
elseif isvector(options.tail)
    options.tail=repmat(options.tail(:)',[1 1 size(options.DES,2)]);
end

if length(var_type)==1
    var_type=var_type.*ones(1,size(Y,2));
end


[Y, orig,w,labels]= expand_categ(Y,var_type,options.labels.dims{2});
options.orig=orig;
    if not(isfield(options,'tail_categ'))|(isfield(options,'tail'))
        n_categ=unique(orig)';
        options.tail_categ=ones(1,0,size(options.tail,3));
        for i=1:length(n_categ)
           options.tail_categ=[options.tail_categ   repmat(options.tail(1,i,:),[1 sum(orig==n_categ(i)) 1])];
        end
    end

options.tail_categ=shiftdim(options.tail_categ,1);

opts=options;
opts.tail=options.tail_categ;

opts.OUT=0;
if ((((comb_funct=='AD'|comb_funct=='A2')|comb_funct=='KS')|comb_funct=='D') | comb_funct=='M')
    opts.Pobs=1;
end
        

opts.Pobs=options.Pobs;
if length(opts.OUT)==1
    opts.OUT=0;
else
    opts.OUT=options.OUT(2);
end

[P_sub T_sub]=NP_ReM(Y,X,DES,B,options.tail_categ,opts.OUT);


%P=ones(size(B,1), length(options.labels.dims{2}),size(DES,1));
%T=zeros(size(B,1), length(options.labels.dims{2}),size(DES,1));
opts.labels.dims{2}=labels(orig==i);
opts.labels.dimslabel{2}=options.labels.dims{2}{1};
if (comb_funct(1)=='A' )|(comb_funct=='KS')
    for i=unique(orig)
        if opts.OUT==1, fprintf(['\n --- Partial p-value for variable: ' options.labels.dims{2}{i} '\n']);end
        if comb_funct(1)=='A'
            opts.w=[];
            for iii=1:size(opts.DES,2)
                ids=zeros(size(X,1),1);
                lbl=unique(X);
                for iv=lbl(find(not(opts.DES(:,iii)==0)))'
                    ids=ids+(X==iv)*sign(opts.DES(iv,iii));
                end
                N=NaNsum(Y(abs(ids)==1,orig==i));
                Npiu=NaNsum(Y(ids==1,orig==i));
                opts.w(1,:,iii)=[Npiu.*(N-Npiu)].^.5;
            end
        elseif comb_funct=='KS'
            opts.w=ones(size(T_sub(1,orig==i,:)));
        end
        if comb_funct=='A2' 
            [P(:,i,:) T(:,i,:)]=NPC(T_sub(:,orig==i,:).^2,'D',opts);
        elseif comb_funct=='KS' 
            [P(:,i,:) T(:,i,:)]=NPC(abs(T_sub(:,orig==i,:)),'M',opts);
        else
            [P(:,i,:) T(:,i,:)]=NPC(T_sub(:,orig==i,:),'D',opts);
        end
    end
elseif comb_funct=='X2',
    
elseif comb_funct=='D' | comb_funct=='M'
    for i=unique(orig)
 %       opts.w=ones(size(T_sub(1,orig==i,i)));
        if opts.OUT==1, fprintf(['\n --- Partial p-value for variable: ' options.labels.dims{2}{i} '\n']);end
        [P(:,i,:) T(:,i,:)]=NPC(T_sub(:,orig==i,:),comb_funct,opts);
    end
else
    for i=unique(orig)
        opts.w=ones(size(P_sub(1,orig==i,:)));
        if opts.OUT==1, fprintf(['\n --- Partial p-value for variable: ' options.labels.dims{2}{i} '\n']);end
        [P(:,i,:) T(:,i,:)]=NPC(P_sub(:,orig==i,:),comb_funct,opts);
    end
end


T=P;
options.p.raw=P(end,:,:,:);
NP_out('Testing repeated measures for categorical variables',options);
