function [P, T, options] = NP_rho(Y,X,Z,B,tail,type,connected,options)
%Nonparametric correlation test for continuous or categorical variables.
%P = NP_rho(Y,X) test the presence of correlation  between any pair of 
%variables, the first in Y and the second in X, where Y is the n-by-m data
%matrix (n sample size, m variables number) of responses and X is another
%n-by-k data matrix (n sample size, M variables number) of responses.
%If X or Y are Struct arrays, the data have a structure as outputted by
%xlsimport.m (see example below). 
%If X or Y are cell arrays, the cells are labels of the variables of X 
%and Y and variables must be declared global variables using reminD.m
%(and have a structure as in xlsimport.m, see example below). 
%P is the (B+1)-by-m-k-C P-values matrix, where C is the number of strata 
%(see below) and B is the number of random permutations set at 1000 by
%default. The element (b,i,j,c) of P refers to b-th permutation of the
%test of correlation between variable i in matrix Y and variable j
%in matrix X of strata c. The last row of P is related to the observed 
%sample.
%
%P = NP_rho(Y,X,Z) perform the correlation test by strata. Z is a vector of
%length n. The C modalities of Z indicate to which strata the
%corresponding rows of Y and X belong. If there are no strata, we can put a 
%vector with a constant value or an empty element (i.e. "[]").
%
%P = NP_rho(Y,X,Z,B) performs B random permutations. If B is a scalar,
%it indicates the number of random permutations. If B is a (B+1)-by-n matrix,
%it indicates the permutation sample space; in this case, each row is a
%(random) permutation of the vector 1:n. The last row must be the vector
%1:n itself if we wish to make inference on the observed data. B = 1000
%by default. 
%
%P = NP_rho(Y,X,Z,B,tail) The default alternative is "correlation is not 0"
%With TAIL it is possible to specify other alternatives:
%
%               TAIL =  0,  alternative: "correlation is not 0".
%               TAIL =  1,  alternative: "correlation is greater than 0".
%               TAIL = -1,  alternative: "correlation is less than 0".
%   
%If TAIL is an m-by-k matrix, the element(i,j) refers to the tail of the
%test of correlation between variable i in matrix Y and variable j
%in matrix X. If TAIL is a scalar, this is considered for each variable.
%
%P = NP_rho(Y,X,Z,B,tail,type) to specify the kind of correlation to test.
%Possible settings for this type are:
%
%       'Pearson'  (the default) to compute Pearson's linear correlation
%                  coefficient (faster if NaNs are not presents). 
%       'Kendall'  to compute Kendall's tau (in this case the Matlab
%                  function corr.m is used)
%       'Spearman' to compute Spearman's rho (in this case the Matlab
%                  function corr.m is used).
%
%The default type of correlation is Pearson’s.
%
%P = NP_rho(Y,X,Z,B,tail,type,connected) to specify the kind of
%connection between strata. Possible settings for the connection are:
%
%               0 to perform independent permutations on strata,
%               1 to perform the same permutations on each stratum.
%
%P = NP_rho(...,'OPTIONS') Possible settings for options are:
%
%	options.OUT =  1 print the observed p-values, 0 don't show it. OUT = 1 by
%                  default.
%   options.Pobs = 1 compute p-values only for observed (not permuted)
%                  data. This saves a considerable amount of time if the number
%                  of variables or permutations is high. Pobs = 0 by default
%                  (compute p-values for both observed and permuted data).
%   options.labels.dims{2} = {'Var1', ......,'Varm'} to customize the
%                  labels of p-values in the output for the m variables Y.
%                  By default the labels are Y1, ..., Ym.
%   options.labels.dims{3} = {'Var1', ......,'Vark'} to customize the
%                  labels of p-values in the output for the q variables X.
%                  By default the labels are X1, ..., Xk.
%   options.labels.dims{4} = {'Strata1', ......,'StrataC'} to customize the
%                  labels of p-values in the output for the C strata.
%                  By default the labels are Strata 1, ..., Strata C.
%
%[P,T] = NP_rho(...)returns the matrix P and the value of the test
%statistics in the (B+1)-by-m-k-C T matrix. The last row is related to the
%observed sample.
%
%[...,OPTIONS] =  NP_rho(...) saves the options used for the analysis. 
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
%The permutation methodology can handle missing data 
%(whether missing at random or not). Please indicate this as NaN.
%
%Example:
%We consider the data set WINE in the matlab file wineRho.mat.
%Below is the data set with four variables: times, Wine, Price,Country. We
%consider the Time variable as a stratification variable.
%
%               Time	Wine  Price	 Country
%                1       3      2      4
%                1       2      3      1
%                1       2      3      1
%                1       3      1      3
%                1       3      1      3
%                1       3      2      4
%                1       3      2      4
%                1       3      2      4
%                1       2      3      1
%                1       3      2      4
%                2       2      3      4
%                2       2      2      2
%                2       2      3      4
%                2       1      2      3
%                2       2      3      4
%                2       2      2      2
%                2       2      3      4
%                2       2      2      2
%                2       1      2      3
%                2       2      2      2
%
%To test the presence of a correlation between the variable Wine and the
%variables Price and Country, strata by strata, using the following 
%vectors as input:
%
%load wineRho
%x = wineRho(:,[3 4]);
%y = wineRho(:,2);
%strata = wineRho(:,1);
%[P,T,options] = NP_rho(y,x,strata);
%
%which produce the output:
%
%load wineRho
%________________________________________________________________
%
% Correlation (mixed moment) test
% p-values: 
%
%_____STRATA    : Strata 1  _____
%			 VARIABLES Y
%VARIABLES X		Y1        
%        X1		0.012987   
%        X2		0.012987   
%
%
%_____STRATA    : Strata 2  _____
%			 VARIABLES Y
%VARIABLES X		Y1        
%        X1		0.46953    
%        X2		1  
%
%To customize the labels of X and Y variables and the strata, we must define
%the variable labels in the options:
%
%options.labels.dims{2} = {'Wine'};
%options.labels.dims{3} = {'Price','Country'};
%options.labels.dims{4} = {'t = 1','t = 2'};
%
%and repeat the test with the options (if options is used, it must
%also specify the B, tail, type and connected):
%
%[P,T,options] = NP_rho(y,x,strata,1000,0,'Pearson',0,options);
%
%which produce the output: 
%
%________________________________________________________________
%
% Correlation (mixed moment) test
% p-values: 
%
%_____STRATA    : t = 1     _____
%			 VARIABLES Y
%VARIABLES X		Wine      
%     Price		0.010989   
%   Country		0.010989   
%
%_____STRATA    : t = 2     _____
%			 VARIABLES Y
%VARIABLES X		Wine      
%     Price		0.46953    
%   Country		1          
%
%These data are also reported in the Excel file wineRho.xls, with column
%labels showing the names of the variables. To perform the same analysis with
%this file, we use the structural arrays D generated with xlsimport:  
%
%D = xlsimport('wineRho');
%[P,T,options] = NP_rho(D(:,2),D(:,[3 4]),D(:,1));
%
%This produces the output:
%
%________________________________________________________________
%
% Correlation (mixed moment) test
% p-values: 
%
%_____STRATA    : Strata 1  _____
%			 VARIABLES Y
%VARIABLES X		Wine      
%     Price		0.005994   
%   Country		0.005994   
%
%_____STRATA    : Strata 2  _____
%			 VARIABLES Y
%VARIABLES X		Wine      
%     Price		0.45255    
%   Country		1          
%
%It is possible to perform the same analysis, with the same output as above,
%using the variable labels of the columns in the Excel files as input.
%Beforehand, it is necessary to declare the variable D global:
%
%D = xlsimport('wineRho');
%reminD(D);
%[P,T,options] = NP_rho({'Wine'},{'Price' 'Country'},'Time');
%
%This function is part of NPClib.
%NPC Lib Matlab library (companion software of F. Pesarin, L. Salmaso: Permutation Tests for Complex Problems: theory, applications and software, John Wiley and sons).
%developed by Livio Finos; consulting team: Rosa Arboretti, Francesco Bertoluzzo, Stefano Bonnini, Chiara Brombin, Livio Corain, Luigi Salmaso and Aldo Solari.
%E-mail for technical questions: livio@stat.unipd.it (Livio Finos).



if nargin <= 1
    error('Not enough iNPut argument')
elseif nargin == 2
   Z=[];
   B=1000;
   tail=zeros(size(Y,2),size(X,2));
   type='Pearson';
   connected=0;
elseif nargin == 3
   B=1000;
   tail=zeros(size(Y,2),size(X,2));
   type='Pearson';
   connected=0;
elseif nargin == 4
   tail=zeros(size(Y,2),size(X,2));
   type='Pearson';
   connected=0;
elseif nargin == 5
   type='Pearson';
   connected=0;
elseif nargin == 6
   connected=0;
end


if nargin==8
    [Y, options]=getDopts(Y,2,options);
else
    [Y, options]=getDopts(Y,2);
end

if isempty(X) X=Y; end
[X, options]=getDopts(X,3,options);

if isempty(Z) Z=ones(size(Y,1),1); end
[Z options]=getDopts(Z,4,options);

if isfield(options,'labels')
 options.labels.dimslabel{4}=options.labels.dims{4};
 temp=[];
 for i=unique(Z)'
    temp=[temp {[options.labels.dims{4}{:} '_' num2str(i)]}];
 end
 options.labels.dims{4}=temp;
end



id=find(not(isnan(Z)));
Z=Z(id);
Y=Y(id,:);
X=X(id,:);
[Z id]=sort(Z);

Y=Y(id,:);
X=X(id,:);
[N NY]=size(Y);
[N NX]=size(X);


N=tabulate(Z);
stratas=N(:,1);
N=N(:,2);
nstratas=length(stratas);

if size(tail,1)>1
    tail=shiftdim(tail,-1);
end
    
if length(tail)==1
    tail=tail.*ones([1 NY NX nstratas]);
else
    if size(tail,2)==1     tail=repmat(tail,[1 NY]); end
    if size(tail,3)==1     tail=repmat(tail,[1 1 NX]); end
    if size(tail,4)<nstratas     tail=repmat(tail,[1 1 1 nstratas]); end
end


if isscalar(B)
    B = Space_perm(N,B,connected);
end 

T=zeros(size(B,1),NY,NX,nstratas);

if not(isempty(strmatch(type,'Spearman','exact'))) &  ((sum(isnan(Y(:)))+sum(isnan(X(:))))==0)
    type='Pearson';
    for ii=1:NX,         X(:,ii)=tiedrank(X(:,ii),0);     end
    for ii=1:NY,         Y(:,ii)=tiedrank(Y(:,ii),0);     end
end

if isempty(strmatch(type,'Pearson','exact'))| (sum(isnan(Y(:)))+sum(isnan(X(:))))
    for j=1:nstratas
        questi=find(Z==stratas(j));
        if length(questi)>1
            Y(:,find(var(Y(questi,:))==0))=NaN;
            X(:,find(var(X(questi,:))==0))=NaN;
            for i=1:size(B,1)
                [RHO]=corr(Y(questi,:),X(B(i,questi),:),'rows','pairwise','type',type);
                T(i,:,:,j)=shiftdim(RHO,-1);
            end
        else
                T(:,:,:,j)=0;
        end
    end
else
    for j=nstratas:-1:1
        questi=find(Z==stratas(j));
        if length(questi)>1
            for i=1:size(B,1)
                T(i,:,:,j)=shiftdim((Y(questi,:)')*X(B(i,questi),:),-1);
            end
            T(:,:,:,j)=T(:,:,:,j)./length(questi);
            mtemp(:,:,j)=(nanmean(Y(questi,:))')*nanmean(X(questi,:));
            %stemp(:,:,j)=(nanstd(Y(questi,:))')*nanstd(X(questi,:));
        else
                T(:,:,:,j)=0;
                mtemp(NY,NX,j)=0;
        end 
    end
    mtemp=shiftdim(mtemp,-1);
    %stemp=shiftdim(stemp,-1);
    for i=1:size(B,1)
        T(i,:,:,:)=(T(i,:,:,:)-mtemp);%./stemp;
    end
end
[P T]=t2p(T,options.Pobs,tail);
options.p.raw=P(end,:,:,:);
NP_out('Correlation (mixed moment) test',options)