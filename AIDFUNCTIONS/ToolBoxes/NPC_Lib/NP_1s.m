function [P, T, options] = NP_1s(Y,B,tail,options)
%Nonparametric permutation test for one sample or paired observations with
%continuous variables.
%P = NP_1s(Y) performs a nonparametric test to compare the means of
%the same sample observed in two different contexts (for example before and
%after treatment). Y is the n-by-m data matrix (n sample size, m variable
%numbers). For a paired sample Y is the differences between the 2 vectors.
%If Y is Struct arrays, the data have a structure as outputted by
%xlsimport.m (see example below). If Y is cell arrays, the cells are labels
%of the variables. The variables must be declared global variables using
%reminD.m (and have a structure as in xlsimport.m, see example below).
%The tested hypothesis tested is "mean=MU"  (MU=0 for each variable by default).
%P is the (B+1)-by-m P-values matrix, where B is the number of random
%permutations set equal to 1000 by default. The last row of P is related to
%the observed sample. The test's default alternative is "mean is not MU"
%The default alternative can be changed by parameter TAIL (see below).
%
%P = NP_1s(Y,B) performs B random permutations. If B is a scalar,
%it indicates the number of random permutations. If B is a (B+1)-by-n matrix,
%it indicates the permutation sample space; in this case, each row is a
%(random) permutation of the vector of +1 and -1. The last row has to be the vector
%1:n itself if we want to make inference on the observed data. B = 1000
%by default. 
%
%P = NP_1s(Y,B,TAIL)is the parameter TAIL and specifies the direction of the 
%alternatives. Possible settings for TAILS are:
%
%	     TAIL =  0, alternative: "mean is not MU".
%        TAIL =  1, alternative: "mean is greater than MU"
%        TAIL = -1, alternative: "mean is less than MU"
%
%TAIL = 0 by default. If TAIL is a 1-by-m vector, each element corresponds
%to a variable; if it is a scalar, this is considered for each variable.
%
%P = NP_1s(...,'OPTIONS') Possible settings for options are:
%
%	options.OUT =  1 print the observed p-values, 0 don't show it. OUT= 1 by
%                  default.
%   options.Pobs = 1 compute p-values only for observed (not permuted)
%                  data. This saves a considerable amount of time if the number
%                  of variables or permutations is high. Pobs = 0 by default
%                  (compute p-values for both observed and permuted).
%   options.labels.dims{2} = {'Var1', ......,'Varm'} to customize the
%                  labels of p-values in the output for the m variables Y.
%                  By default the labels are Y1, ..., Ym.
%   options.MU     mean considered for the testing problem. If MU is a 1-by-m
%                  vector, each element corresponds to a variable; if it is a
%                  scalar, this is considered for each variable. MU = 0 by
%                  default.
%
%[P,T] = NP_1s(...) returns the P matrix and value of the test statistics
%in the (B+1)-by-m T matrix. The last row is related to the observed
%sample.
%
%[...,OPTIONS] =  NP_1s(...) saves the options used for the analysis. 
%OPTIONS is a structural array with the following structure:
%   
%       options.labels.dimslabel label for the dimensions of the P-value
%                matrix P:
%		         First dimension label 'Random Permutation'
%		         Second dimension label 'VARIABLES Y'
%	    options.labels.dims labels for any variable of the columns of the
%                P-value matrix P
%       options.OUT see above
%       options.Pobs see above
%	    options.tail see above
%       options.Combdims dimension of the combination. By default it is set
%                equal to length(size(P)), i.e. the last dimension of 
%                matrix P
%	    options.p.raw    P-values
%
%Example
%As first example we consider the data set blood1s in the matlab file
%blod1s.mat. Below is the data set of the differences of a univariate blood
%variable observed three times after two different treatments in ten
%patients. 
%
%           Patient     t = 1	t = 2	t = 3
%               1        0,1	-0,4	-0,8
%               2        0,1	-0,4	-0,6
%               3       -0,2	-0,5	-0,1
%               4       -0,4	-0,8	-0,4
%               5        0,2	 0,1     0,2
%               6       -0,2	-1,1	-0,4
%               7        0,3	-0,3	-0,7
%               8       -0,3	 0,3	 0,3
%               9        0       0,2	-0,4
%               10      -0,5	-0,7	-0,9
%
%To perform a one-sample test with paired observations on these data, the
%following vectors are used as input:
%
%load blood1s;
%y = blood1s(:,2:4);
%[P,t,options] = NP_1s(y);
%
%which produce the output:
%
%________________________________________________________________
%
% Testing Zero-mean distribution for 1 sample 
% for continuous (or dichotomous) variables
% p-values: 
%			        VARIABLES Y
%           		Y1        	Y2        	Y3        
%           		0.3956     	0.047952   	0.027972
%
%To customize the labels of Y variables, we must define
%the variable labels in the options as in the previous example:
%
%options.labels.dims{2} = {'t = 1','t = 2','t = 3'};
%
%and repeat the test with the options (if options is used, the B and TAIL
%parameters must also be specified):
%
%[P,t,options] = NP_1s(y,1000,0,options);
%
%which produce the output: 
%
%________________________________________________________________
%
% Testing Zero-mean distribution for 1 sample 
% for continuous (or dichotomous) variables
% p-values: 
%                   VARIABLES Y
%           		t = 1     	t = 2     	t = 3     
%           		0.37662    	0.036963   	0.025974  
%
%These data are also reported in the Excel file blood1s.xls, with column
%labels showing the name of the variables. To perform the same analysis with this
%file we use the structural arrays D generated with xlsimport: 
%
%[D]=xlsimport('blood1s');
%[P,T,options] = NP_1s(D(:,2:end));
%
%This produces the output:
%
%________________________________________________________________
%
% Testing Zero-mean distribution for 1 sample 
% for continuous (or dichotomous) variables
% p-values: 
%                   VARIABLES Y
%           		t = 1     	t = 2     	t = 3     
%           		0.38162    	0.042957   	0.023976 
%
%It is possible to perform the same analysis, with the same output as above,
%using the variable labels of the columns in the Excel files as input.
%Beforehand, the variable D must be declared global:
%
%[D]=xlsimport('blood1s');
%reminD(D)
%[P,T,options] = NP_1s({'t = 1','t = 2','t = 3'});
%
%This function is part of NPClib.
%NPC Lib Matlab library (companion software of F. Pesarin, L. Salmaso: Permutation Tests for Complex Problems: theory, applications and software, John Wiley and sons).
%developed by Livio Finos; consulting team: Rosa Arboretti, Francesco Bertoluzzo, Stefano Bonnini, Chiara Brombin, Livio Corain, Luigi Salmaso and Aldo Solari.
%E-mail for technical questions: livio@stat.unipd.it (Livio Finos).

if nargin==4
    [Y, options]=getDopts(Y,2,options);
else
    [Y, options]=getDopts(Y,2);
end
    options=get_options(ones(1,size(Y,2),size(Y,3)),'1s',options);

if nargin == 1    
   tail=zeros(1,size(Y,2));
   B=1000;
elseif nargin == 2
   B=1000;
end

if length(tail)==1
    tail=tail.*ones(1,size(Y,2));
end


if length(options.MU)==1
    options.MU=options.MU.*ones(1,size(Y,2));
end


[n m]=size(Y);

Y=Y-options.MU(ones(n,1),:);

if prod(size(B))==1 & B>0
    
T((B+1),:)=  nanmean(Y);%./std(Y) ;
for i=1:B
   perm=repmat(2.*binornd(1,.5,n,1)-1,1,m);
   T(i,:)= nanmean(perm.*Y);%./std(perm.*Y) ;
end

else
    if prod(size(B))==1  %cioè B==0
        B=space_sign(n,0);
    end
      perm=B';
    BB=size(B,1);
    B=BB-1;
    T=zeros(BB,m);
for i=1:BB
   T(i,:)= nanmean(repmat(perm(:,i),1,m).*Y);%./std(repmat(perm(:,i),1,m).*Y) ;
end
end

[P T]=t2p(T,options.Pobs,tail);

        options.p.raw=P(end,:,:);
    NP_out('Testing Zero-mean distribution for 1 sample \n for continuous (or dichotomous) variables\n p-values:',options)
