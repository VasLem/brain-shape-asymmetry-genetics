function [P, Chi2, options, Partial]  = NP_Chi2(Y,X,B,options)
%Nonparametric independence test for categorical variables.
%P = NP_Chi2(Y,X) test the independence using the Chi2 statistic on the
%conjoint frequencies table between any pair of categorical variables, the
%first in Y and the second in X, where Y is the n-by-m data matrix (n sample
%size, m variables number) of responses and X is another n-by-M data matrix
%(n sample size, M variables number) of responses.
%If X or Y are Struct arrays, the data have a structure as outputted by
%xlsimport.m. 
%If X or Y are cell arrays, the cells are labels of the variables of X 
%and Y and variables must be declared global variables using reminD.m
%(and have a structure as in xlsimport.m). 
%
%P = NP_Chi2(Y,X,B)  performs B random permutations. If B is a scalar,
%it indicates the number of random permutations. If B is a (B+1)-by-n matrix,
%it indicates the permutation sample space; in this case, each row is a
%(random) permutation of the vector 1:n. The last row have to be the vector
%1:n itself if you want to make inference on the observed data. B = 1000
%by default. 
%
%P = NP_Chi2(...'options') Possible settings for options are:
%
%	options.OUT =  1 print the observed p-values, 0 don't show it. OUT = 1 by
%                  default.
%   options.Pobs = 1 compute p-values only for observed (not permuted)
%                  data. This save a considerable amount of time if the number
%                  of variables or permutations is high. Pobs = 0 by default
%                  (compute p-values for both observed and permutated data).
%   options.labels.dims{2} = {'Var1', ......,'Varm'} to customize the
%                  labels of p-values in the output for the m variables Y.
%                  By defoult the labels are Y1, ..., Ym
%   options.labels.dims{3} = {'Var1', ......,'Varq'} to customize the
%                  labels of p-values in the output for the q variables X.
%                  By default the labels are X1, ..., Xq
%
%[P Chi2] = NP_Chi2(...) also returns the (B+1)-by-1 Chi2 statistics
%vector. The last row is related to the observed sample.
%
%[...,options] NP_Chi2(...) saves the options used for the analysis. 
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
%This function is part of NPClib.
%NPC Lib Matlab library (companion software of F. Pesarin, L. Salmaso: Permutation Tests for Complex Problems: theory, applications and software, John Wiley and sons).
%developed by Livio Finos; consulting team: Rosa Arboretti, Francesco Bertoluzzo, Stefano Bonnini, Chiara Brombin, Livio Corain, Luigi Salmaso and Aldo Solari.
%E-mail for technical questions: livio@stat.unipd.it (Livio Finos).


if nargin <= 1
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


    

[N m]=size(Y);
if isscalar(B)
    [B] = Space_perm(N,B);
end


x=X;
y=Y;
for xi=size(X,2):-1:1
X=x(:,xi);

X= expand_categ(X(:),1);
[Y, orig]= expand_categ(y,1);
ids=unique(orig);

expect=sum(Y)'*sum(X)./N;

    perm=B';
    
    for b=size(B,1):-1:1
        
        Partial(b,:,:,xi)=(Y(perm(:,b),:)'*X-expect).^2./expect;
    end
    for i=1:m
        Chi2(:,i,xi)=sum(sum(Partial(:,orig==ids(i),:),3),2);
    end
    
end
    
    [P Chi2]=t2p(Chi2,options.Pobs);

if options.OUT==1
    options.p.raw=reshape(P(end,:,:),size(P(end,:,:),2),size(P(end,:,:),3))';
    NP_out('Testing independence of Ys with respect to X \n for chategorical variables \n p-values',options)
     fprintf(' Chi2:         ')
     fprintf('\t%4.4f  ',Chi2(end,:))
     fprintf('\n')
end