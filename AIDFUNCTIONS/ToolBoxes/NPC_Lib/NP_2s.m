function [P, T, options] = NP_2s(Y,X,B,tail,options)
%Nonparametric permutation one-way ANOVA test for equality in distribution
%of 2 samples.
%P = NP_2S(Y,X) performs a one-way ANOVA to compare the means of two 
%samples. Y is the n-by-m data matrix (n=n1+n2 sample size, 
%m variables number) of responses. X is a vector of length n. The two 
%modalities of X indicate to which sample the corrisponding row of Y belongs.
%If X is an n-by-q matrix, the same analysis is replicated for each column.
%If X or Y are struct arrays, the data have a structure as outputed by
%xlsimport.m (see example below). 
%If X or Y are cell arrays, the cells are labels of the variables of X 
%and Y and variables must be declared global variables using reminD.m
%(and have a structure as in xlsimport.m, see example below). 
%P is the (B+1)-by-m-by-q P-values matrix, where B is the number of
%random permutations set equal to 1000 by default. The last row of P is
%related to the observed sample. The default alternative of the test is 
%"means are not equal". The default alternative can be changed by
%parameter TAIL (see below)
%
%P = NP_2S(Y,X,B) performs B random permutations. If B is a scalar,
%it indicates the number of random permutations. If B is a (B+1)-by-n matrix,
%it indicates the permutation sample space; in this case, each row is a
%(random) permutation of the vector 1:n. The last row have has be the vector
%1:n itself if we wish to make inference on the observed data. B = 1000
%by default. 
%
%P = NP_2S(Y,X,B,TAIL), the default alternative is "means are not equal". 
%With TAIL it is possible to specify other alternatives: 
%
%	     TAIL =  0, alternative: "means are not equal"
%        TAIL =  1, alternative: "mean of first group is less than mean
%                   of the second group"
%        TAIL = -1, alternative: "mean of first group is greater than
%                   mean of the second group"
%
%If TAIL is a vector (or a matrix if X is a matrix), each element
%corresponds to a variable; if it is  a scalar, this is considered for each
%variable.
%
%P = NP_2S(...,'OPTIONS') Possible settings for options are:
%
%	options.OUT =  1 print the observed p-values, 0 don't show it. OUT= 1 by
%                  default.
%   options.Pobs = 1 compute p-values for observed (not permutated) data
%                  only, it saves a considerable amount of time if the number
%                  of variables or permutations is high. Pobs = 0 by default
%                  (compute p-values for observet data as well as for
%                  permuted).
%   options.labels.dims{2} = {'Var1', ......,'Varm'} to customize the
%                  labels of p-values in the output for the m variables Y.
%                  By default the labels are Y1, ..., Ym
%   options.labels.dims{3} = {'Var1', ......,'Varq'} to customize the
%                  labels of p-values in the output for the q variables X.
%                  By default the labels are X1, ..., Xq
%
%[P,T] = NP_2S(...) returns the P matrix and value of the test statistics
%in the (B+1)-by-m-by-q T matrix. The last row is related to the observed
%sample.
%
%[...,OPTIONS] =  NP_2S(...) saves the options used for the analysis. 
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
%	    options.tail see above
%       options.Combdims dimension of the combination. By default it is set 
%                equal to length(size(P)), i.e. the last dimension of
%                matrix P
%	    options.p.raw    P-values
%
%The permutation methodology can handle missing data 
%(whether missing at random or not). Please indicate this as NaN.
%
%Example: As a first example we consider the data set BLOODPR in the matlab
%file bloodpr.mat. Below is the data set with two samples (GROUP 1 and GROUP 2)
%and one variable (RESPONSE):
%
%               GROUP	RESPONSE
%               1       94
%               1       108
%               1       110
%               1       90
%               2       80
%               2       94
%               2       85
%               2       90
%               2       90
%               2       90
%               2       108
%               2       94
%               2       78
%               2       105
%               2       88
%
%To perform a two-sample test with these data, the following vectors are
%used as input:
%
%load bloodpr 
%x = bloodpr(:,1);
%y = bloodpr(:,2);
%[P,T,options] = NP_2s(y,x);
%
%which produce the output:
%
%                   VARIABLES Y
%           		Y1  
%                   0.0999 
%
%To customize the labels of X and Y variables, we must define
%the variable labels in the options:
%
%options.labels.dims{2} = {'RESPONSE'};
%options.labels.dims{3} = {'GROUP'};
%
%and repeat the test with the options (if options is used, the B and TAIL
%parameters must also be specified):
%
%[P,T,options] = NP_2s(y,x,B,TAIL,options);
%
%which produce the output: 
%
%                   VARIABLES Y
%           	    RESPONSE  
% GROUP		        0.0999   
%
%We consider the more complicated two-sample, four variables example
%"Germina" discussed in the reference book (Note: there are also missing
%data in this data set). The data are saved in a binary "MAT-file" 
%germina.mat. To perform a two-sample test between the Fertlized sample 
%and the Non Fertilized sample with the variables Germinated, Weight,
%Surface and Surface2:  
%
%load germina 
%x = germina(:,1);
%y = germina(:,2:end);
%[P,T,options] = NP_2s(y,x); 
%
%which produce the output:
%
%                   VARIABLES Y
%           		Y1        	Y2        	Y3        	Y4        
%           		0.17083    	0.17982    	0.041958   	0.065934 
%
%To customize the labels of X and Y variables, we must define
%the variable labels in the options as in the previous example:
%
%options.labels.dims{2} = {'Germinate','Weight','Surface','Surface2'};
%options.labels.dims{3} = {'Fertilizer'};
%[P,T,options] = NP_2s(y,x,B,TAIL,options); 
%
%These data are also reported in the Excel file germina.xls, with column
%labels showing the names of the variables. To perform the same analysis
%whit this file we use the structural arrays D generated with xlsimport:  
%
%[D]=xlsimport('germina');
%[P,T,options] = NP_2s(D(:,2:5),D(:,1));
%
%This produces the output:
%
%               VARIABLES Y
%           	Germinated	Weight    	Surface   	Surface2  
%Fertilizer		0.17083    	0.17982    	0.041958   	0.065934
%
%It is possible to perform the same analysis, with the same output as above,
%using the variable labels of the columns in the Excel files as input.
%Beforehand, the variable D must be declared global:
%
%[D]=xlsimport('germina');
%reminD(D)
%[P,T,options] = NP_2s({'Germinated','Weight','Surface','Surface2'},'Fertilizer');
%
%
%This function is part of NPClib.
%NPC Lib Matlab library (companion software of F. Pesarin, L. Salmaso: Permutation Tests for Complex Problems: theory, applications and software, John Wiley and sons).
%developed by Livio Finos; consulting team: Rosa Arboretti, Francesco Bertoluzzo, Stefano Bonnini, Chiara Brombin, Livio Corain, Luigi Salmaso and Aldo Solari.
%E-mail for technical questions: livio@stat.unipd.it (Livio Finos).


if nargin < 2
    error('Not enough iNPut argument')
elseif nargin == 2
   B=1000;
   tail=0;
elseif nargin == 3
   tail=0;
end

if nargin==5
    [Y, options]=getDopts(Y,2,options);
else
    [Y, options]=getDopts(Y,2);
end

[X, options]=getDopts(X,3,options);
    sizY=size(Y,2);
    sizX=size(X,2);
    options=get_options(ones([1 sizY sizX]),'partial',options);
options.tail=tail;

[origsizeY]=size(Y);
N=origsizeY(1);
m=prod(origsizeY(2:end));
Y=reshape(Y,N,m);

    cc=unique(X)';
NOTmiss=isfinite(Y);

if isscalar(B)
    [B] = Space_perm(N,B);
end

T=zeros(size(B,1),m,size(X,2));
options.B=size(B,1);

Y(find(1-NOTmiss))=0;
   % sum_Y=sum(Y);
    
    mancanti=find(sum(NOTmiss)<size(Y,1));
    nonmanc=setdiff((1:size(Y,2)),mancanti);
    
for j=1:size(X,2)
    x=X(:,j);
    cc=unique(x)';
    if length(cc)>2
        fprintf('\nWARNING: NP_2s works only with 2 samples, but you are iNPuting %1.0f samples. \n Use NP_ANOVA1 instead.\n',length(cc))
    elseif length(cc)==1
            T(1:size(B,1),1:size(Y,2),j)=0;
    else
    
    n=[(sum(x==cc(1))) (sum(x==cc(2)))];
    
    nu=sum(NOTmiss(:,mancanti));
    
        BB=size(B,1);
        if not(isempty(nonmanc))
            deno(find(x==cc(2)))=1/sum(x==cc(2));
            deno(find(x==cc(1)))=-1/sum(x==cc(1));
            deno=deno(:)';
        for i=1:BB
            T(i,nonmanc,j)=deno*Y(B(i,:),nonmanc);
        end
        end
        if not(isempty(mancanti))
        for i=1:BB
            nu_samp=sum(NOTmiss(B(i,x==cc(1)),mancanti),1);
            T(i,mancanti,j)= sum(Y(B(i,x==cc(2)),mancanti),1).*(    nu_samp ./(nu-nu_samp+(nu_samp==nu) )).^.5-...
                sum(Y(B(i,x==cc(1)),mancanti),1) .*((nu-nu_samp)./(nu_samp   +(nu_samp==0)  )).^.5;
        end
        end
        
        end
end
T=reshape(T,[size(T,1) origsizeY(2:end) size(T,3)]);
[P T]=t2p(T,options.Pobs,options.tail);

        options.p.raw=P(end,:,:);
 NP_out('Testing equality in distribution for 2 independent samples \n for continuous (or dichotomous) variables (with or without missing values)',options)
