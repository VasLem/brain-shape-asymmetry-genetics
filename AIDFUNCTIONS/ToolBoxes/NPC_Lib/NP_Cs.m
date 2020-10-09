function [P, T, options] = NP_Cs(Y,X,B,options)
%NONPARAMETRIC permutation test for C > 2 samples for continuous variables.
%P = NP_Cs(Y,X) performs a one-way ANOVA to compare the means of C > 2 
%samples. Y is the n-by-m data matrix (n=n1+n2+ ... +nC sample size, 
%m variables number) of responses. X is a vector of length n. The C 
%modalities of X indicate to which sample the corrisponding row of Y belongs.
%If X is an n-by-q matrix, the same analysis is replicated for each column.
%If X or Y are Struct arrays, the data have a structure as outputted by
%xlsimport.m (see example below). 
%If X or Y are cell arrays, the cells are labels of the variables of X 
%and Y and variables must be declared global variables using reminD.m
%(and have a structure as in xlsimport.m, see example below). 
%P is the (B+1)-by-m-by-q P-values matrix, where B is the number of
%random permutations set equal to 1000 by default. The last row of P is
%related to the observed sample. The alternative of the test is "means are not
%equal". 
%
%P = NP_Cs(Y,X,B) performs B random permutations. If B is a scalar,
%it indicates the number of random permutations. If B is a (B+1)-by-n matrix,
%it indicates the permutation sample space; in this case, each row is a
%(random) permutation of the vector 1:n. The last row must be the vector
%1:n itself if we wish to make inference on the observed data. B = 1000
%by default. 
%
%P = NP_Cs(...,'OPTIONS') Possible settings for options are:
%
%	options.OUT =  1 print the observed p-values, 0 don't show it. OUT = 1 by
%                  default.
%   options.Pobs = 1 compute p-values only for observed (not permuted)
%                  data. This save a considerable amount of time if the number
%                  of variables or permutations is high. Pobs = 0 by default
%                  (compute p-values for for both observed and permuted data).
%   options.labels.dims{2} = {'Var1', ......,'Varm'} to customize the
%                  labels of p-values in the output for the m variables Y.
%                  By default the labels are Y1, ..., Ym
%   options.labels.dims{3} = {'Var1', ......,'Varq'} to customize the
%                  labels of p-values in the output for the q variables X.
%                  By default the labels are X1, ..., Xq
%
%[P,T] = NP_Cs(...) returns the P matrix and value of the test statistics
%in the (B+1)-by-m-by-q T matrix. The last row is related to the observed
%sample.
%
%[...,OPTIONS] =  NP_Cs(...) saves the options used for the analysis. 
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
%The permutation methodology can handle with missing data 
%(whether missing at random or not). Please indicate this as NaN.
% 
%Example:
%We consider the data set PERCRAWDATA in the matlab file percrawdataCs.mat.
%Below the data set with four samples (identified by the variable dose which
%assumes values 0, 150, 1500, 5000) and four response variables (foot_splay,
%weight, temperature and motor_activity). The total number of subject is 25:
%
%   rat     dose	foot_splay	weight	temperature	motor_activity
%   1       150     20.5        140.7       37.6        34
%   2       1500	35          144.6       38.1        13
%   3       150     25.5        151.7       37.8        174
%   4       150     22          161.1       38          210
%   5       0       31          144.8       37.3        64
%   6       5000	43          143.8       34.2        46
%   7       0       32          154         37.5        220
%   8       150     34          152.1       37.6        21
%   9       1500	26          147         37.6        146
%   10      500     39.5        161.2       37.8        96
%   11      1500	39.5        147.6       37.2        137
%   12      0       21          153.4       37.6        81
%   13      5000    NaN         147.2       26.8        0
%   14      150     33          151.5       37.2        216
%   15      1500	35.5        149.9       37.5        154
%   16      150     39          160.4       37.4        221
%   17      1500	36          139.2       37.3        45
%   18      0       34          148.6       37.5        146
%   19      0       29          158.8       37.8        82
%   20      0       32.5        151         37.9        155
%   21      500     36.5        144.5       38          35
%   22      150     36.5        150.9       38.4        79
%   23      500     34          142.2       37.5        33
%   24      500     38.5        151         36.8        168
%   25      1500	27.5        152         38.1        74
%
%To perform a C-sample test with these data, the following vectors are used
%as input:
%
%load percrawdataCs; 
%x = percrawdataCs(:,2);
%y = percrawdataCs(:,3:end);
%[P,T,options] = NP_Cs(y,x);
%
%which produce the output:
%
%                   VARIABLES Y
%           		Y1        	Y2        	Y3        	Y4        
%           		0.12887    	0.34865    	0.002997   	0.3027  
%
%To customize the labels of X and Y variables, we must define
%the variable labels in the options:
%
%options.labels.dims{2} = {'foot_splay','weight','temperature', ...
%        'motor_activity'};
%options.labels.dims{3} = {'dose'};
%
%and repeat the test with the options (if options is used, it must
%also specify the B and TAIL parameters):
%
%[P,T,options] = NP_Cs(y,x,B,options);
%
%which produce the output: 
%
%                   VARIABLES Y
%           		foot_splay	weight    	temperature	moto..ivity
%        dose   	0.12887    	0.34865    	0.002997   	0.3027
%
%These data are also reported in the Excel file percrawdataCs.xls, with
%column labels showing the names of the variables. To perform the same
%analysis with this file we use the structural arrays D generated
%with xlsimport:  
%
%[D]=xlsimport('percrawdataCs');
%[P,T,options] = NP_Cs(D(:,3:end),D(:,2));
%
%This produce the output:
%                   VARIABLES Y
%           		foot_splay	weight    	temperature	moto..ivity
%        dose   	0.12887    	0.34865    	0.002997   	0.3027
%
%It is possible to perform the same analysis, whith the same output as above,
%using the variable labels of the columns in of the Excel files as input.
%Beforehand, it is necessary to declare the variable D global:
%
%[D]=xlsimport('percrawdataCs');
%reminD(D)
%[P,T,options] = NP_cs({'foot_splay','weight','temperature', ...
%   'motor_activity'},'dose');
%
%This function is part of NPClib.
%NPC Lib Matlab library (companion software of F. Pesarin, L. Salmaso: Permutation Tests for Complex Problems: theory, applications and software, John Wiley and sons).
%developed by Livio Finos; consulting team: Rosa Arboretti, Francesco Bertoluzzo, Stefano Bonnini, Chiara Brombin, Livio Corain, Luigi Salmaso and Aldo Solari.
%E-mail for technical questions: livio@stat.unipd.it (Livio Finos).


if nargin < 2
    error('Not enough iNPut argument')
elseif nargin == 2
   B=1000;
end

if nargin==4
    [Y, options]=getDopts(Y,2,options);
else
    [Y, options]=getDopts(Y,2);
end
[X, options]=getDopts(X,3,options);

sizY=size(Y,2);
sizX=size(X,2);
options=get_options(ones([1 sizY sizX]),'partial',options);



[origsizeY]=size(Y);
N=origsizeY(1);
m=prod(origsizeY(2:end));
Y=reshape(Y,N,m);  
miss=isfinite(Y);

if isscalar(B)
    [B] = Space_perm(N,B);
end
T=zeros(size(B,1),m,size(X,2));
options.B=size(B,1);    
Y(find(1-miss))=0;
    sum_Y=sum(Y);
    
    mancanti=find(sum(miss)<size(Y,1));
    nonmanc=setdiff((1:size(Y,2)),mancanti);
    %Ymiss=Y(:,mancanti);
    
  %Y(:,nonmanc)=mncn(Y(:,nonmanc));
for j=1:size(X,2)
    x=X(:,j);
    cc=unique(x)';
    n=zeros(1,0);
    for c=cc
        n=[n (sum(x==c))];
    end

        perm=B';
        BB=size(B,1);
        if not(isempty(nonmanc))
        
        for i=1:BB
            for c=1:length(cc)
                %T(i,nonmanc,j)=T(i,nonmanc,j)+var(Y(perm(x==cc(c),i),nonmanc)).*n(c);
                %T(i,nonmanc,j)=T(i,nonmanc,j)+(sum(Y(perm(x==cc(c),i),nonmanc).^2)-sum(Y(perm(x==cc(c),i),nonmanc)).^2./n(c));
                T(i,nonmanc,j)=T(i,nonmanc,j)+(sum(Y(perm(x==cc(c),i),nonmanc)).^2./n(c));
            end
        end
        
        %T(:,nonmanc,j)=-T(:,nonmanc,j);
        end
        if not(isempty(mancanti))
    nu=sum(miss(:,mancanti));
        for i=1:BB
            for c=1:length(cc)
                nu_samp=sum(miss(perm(x==cc(c),i),mancanti),1);
                T(i,mancanti,j)=T(i,mancanti,j)+ (  sum(Y(perm(x==cc(c),i),mancanti),1) .*((nu-nu_samp)./(nu_samp+(nu_samp==0) )).^.5  -...
                                   (sum_Y(1,mancanti)-sum(Y(perm(x==cc(c),i),mancanti),1)).*(   nu_samp./(nu-nu_samp+(nu_samp==nu) )).^.5 ).^2;
            end
             
        end
        end
%     T(:,nonmanc,j)=repmat((sum(Y(perm(:,i),nonmanc).^2)-sum(Y(perm(:,i),nonmanc)).^2./(N-1)),size(T,1),1)-T(:,nonmanc,j);

end

%T=ceil(T.*10.^ceil(-log10(T)));
T=reshape(T,[size(T,1) origsizeY(2:end) size(T,3)]);
[P T]=t2p(T,options.Pobs);

    options.p.raw=P(end,:,:);
    NP_out('Testing equality in distribution for C independent samples \n for continuous (or binomial) variables (with or without missing values)',options)
