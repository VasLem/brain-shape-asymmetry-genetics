function [P, T, options] = NP_Cs_MA(Y,X,B,stats,comb_funct,options)
%Nonparametric permutation MULTI-ASPECT one-way ANOVA test for equality in
%distribution of C > 2 samples.
%P = NP_Cs_MA(Y,X) performs a one-way ANOVA to compare the means of C > 2 
%samples. Y is the n-by-m data matrix (n=n1+n2+ ... +nC sample size, 
%m variables number) of responses. X is a vector of length n. The C 
%modalities of X indicate to which sample the corrisponding row of Y belongs.
%If X is an n-by-q matrix, the same analysis is replicated for each column.
%If X or Y are Struct arrays, the data have a structure as outputted by
%xlsimport.m.
%If X or Y are cell arrays, the cells are labels of the variables of X 
%and Y and variables must be declared global variables using reminD.m
%(and have a structure as in xlsimport.m). 
%P is the (B+1)-by-m-by-q P-values matrix, where B is the number of
%random permutations set at 1000 by default. The last row of P is
%related to the observed sample. The alternative of the test is "means are not
%equal". 
%
%P = NP_Cs_MA(Y,X,B) performs B random permutations. If B is a scalar,
%it indicates the number of random permutations. If B is a (B+1)-by-n matrix,
%it indicates the permutation sample space; in this case, each row is a
%(random) permutation of the vector 1:n. The last row must be a vector
%1:n itself if we wish to make inference on the observed data. B = 1000
%by default. 
%
%P = NP_Cs_MA(Y,X,B,stats) to perform a MULTI-ASPECT one-way ANOVA test. stats
%is a 1-by-s cell array, where s is the number of considered aspects. You
%Possible settings for stats are:
%
%       'F'                         for standard ANOVA statistics
%       'Kruskal-Wallis' (or 'KW')  for rank-based statistics
%       'Mean_grps'      (or 'M' )  for means comparisons
%       'Median_grps'    (or 'Me')  for medians comparisons
%       'Kolmog-Smirn'   (or 'KS')  for Kolmogorov-Smirnov statistics
%       'Anders-Darl'    (or 'AD'}  for Anderson-Darling statistics.
%
%We can also perform direct transformation on the Ys using a ':'. For
%example:  
%
%   stats={':Y.^2',':R(Y).^2','MW','KS'}
%
%performs the following test:
%
%               1) over the Y's second moments;
%               2) over the ranks;
%               3) using F statistic; 
%               4) using Kruskal-Wallis and Kolmogorov-Smirnov statistics.
%
%P = NP_Cs_MA(Y,X,B,stats,comb_funct) to specify the combining
%function used to combine the test statistics or p-values into a global
%p-value. See NPC for details.
%
%P = NP_2s_MA(Y,X,B,stats,comb_funct,tail) The default alternative is
%"means are not equal". With TAIL it is possible to specify other alternatives: 
%
%	     TAIL =  0, alternative: "means are not equal"
%        TAIL =  1, alternative: "mean of first group is less than mean
%                   of the second group"
%        TAIL = -1, alternative: "mean of first group is greater than
%                   mean of the second group"
%
%If TAIL is a vector (or a matrix if X is a matrix), each element
%corresponds to a variable; if it is a scalar, this is considered for each
%variable.
%
%P = NP_2s_MA(...,'options') Possible settings for options are:
%
%	options.OUT =  1 print the observed p-values, 0 don't show it. OUT= 1 by
%                  default.
%   options.Pobs = 1 compute p-values only for observed (not permutted)
%                  data. This save a considerable amount of time if the number
%                  of variables or permutations is high. Pobs = 0 by default
%                  (compute p-values for both observed and permuted data).
%   options.labels.dims{2} = {'Var1', ......,'Varm'} to customize the
%                  labels of p-values in the output for the m variables Y.
%                  By default the labels are Y1, ..., Ym
%   options.labels.dims{3} = {'Var1', ......,'Varq'} to customize the
%                  labels of p-values in the output for the q variables X.
%                  By default the labels are X1, ..., Xq
%
%
%[P T] = NP_2s_MA(...) returns the P matrix and value of the test statistics
%in the (B+1)-by-m-by-q T matrix. The last row is related to the observed
%sample.
%
%[...,OPTIONS] = NP_2s_MA(...) saves the options used for the analysis. 
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
%The permutation methodology can handle missing data (whether missing at
%random or not). Please indicate this as NaN.
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
if nargin==6
    [Y, options]=getDopts(Y,2,options);
else
    [Y, options]=getDopts(Y,2);
end
[X, options]=getDopts(X,3,options);

    sizY=size(Y,2);
    sizX=size(X,2);
    options=get_options(ones([1 sizY sizX]),'partial',options);


    [N m]=size(Y);
    
miss=isfinite(Y);

if isscalar(B)
    [B] = Space_perm(N,B);
end
    mancanti=find(sum(miss)<size(Y,1));
    nonmanc=setdiff((1:size(Y,2)),mancanti);
    
    opts.OUT=0;
if not(isempty(mancanti))
    [Pm Tm]=NP_Cs(Y(:,mancanti),X,B,opts);
end

Y=Y(:,nonmanc);
[N m1]=size(Y);
P=zeros(size(B,1),m1,length(stats)+1);
T=P;
for j=1:length(stats)
    switch stats{j}(1)
        case {':'}
            stats{j}=stats{j}(2:end);
%            ['[P(:,nonmanc,j) T(:,nonmanc,j)]=NP_Cs(' stats{j} ',X,B,opts);']
            eval(['[P(:,:,j) T(:,:,j)]=NP_Cs(' stats{j} ',X,B,opts);']);
        otherwise
            [P(:,:,j) T(:,:,j)]=NP_Cs_stat(Y,X,B,stats{j});
    end
end

 for j=1:m1
     [P(:,j,end) T(:,j,end)]=NPC(squeeze(P(:,j,1:end-1)),comb_funct,opts);
 end
 stats{end+1}='Global';

 if not(isempty(mancanti))
     stats{end+1}='Miss.Vals.';
     p=ones(size(P,1),m,length(stats)).*NaN;
     p(:,mancanti,end)=Pm;
     p(:,nonmanc,1:end-1)=P;
     P=p;
 end
options.p.raw=P(end,:,:);
options.labels.dimslabel{3}='STATISTIC  ';
options.labels.dims{3}=stats;
    NP_out('Testing equality in distribution for C independent samples \n for continuous (or dichotomous) variables (with or without missing values)\n p-values:',options)



%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y=R(Y)

for i=1:size(Y,2)
    [Y(:,i)] = tiedrank(Y(:,i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [P, T]=NP_Cs_stat(Y,X,perm,stat)
cc=unique(X)';
n=zeros(1,0);
    for c=cc
        n=[n (sum(X==c))];
    end

switch stat
    case {'Median_grps','Me'}
        for i=size(perm,1):-1:1
            for c=1:length(cc)
                medians(c,:)=median(Y(perm(i,X==cc(c)),:)).^2.*n(c);
            end
            T(i,:)=sum(medians);
        end
        [P T]=t2p(T);
    case {'Mean_grps','M'}
        for i=size(perm,1):-1:1
            for c=1:length(cc)
                means(c,:)=mean(Y(perm(i,X==cc(c)),:)).^2.*n(c);
            end
            T(i,:)=sum(means);
        end
        [P T]=t2p(T);
    case {'Kruskal-Wallis','KW','F'}
        if (strcmp(stat,'Kruskal-Wallis')| strcmp(stat,'KW')),Y=R(Y);end
        displayopt='off';
        for i=size(perm,1):-1:1
            for j=size(Y,2):-1:1
                T(i,j)=anova1(Y(:,j),X(perm(i,:)),displayopt);
            end
        end
        T=finv(1-T,length(n),length(X)-length(n)-1);
        [P T]=t2p(T);
    case {'Kolmog-Smirn','KS'}
        [YY I]=sort(Y);
        temp=((1:length(X)))./length(X);
        for i=size(Y,2):-1:1
            t_ecdf(I(:,i),i)=temp;
        end
        
        for i=size(perm,1):-1:1
            p_ecdf=NaN.*t_ecdf;
            for c=cc
                chi=find(X(perm(i,:))==c);
                [YY I]=sort(Y(chi,:));
                temp=(1:length(chi))./length(chi);
                for ii=size(Y,2):-1:1
                    p_ecdf(chi(I(:,ii)),ii)=temp;
                end
            end
            T(i,:)=max(abs(p_ecdf-t_ecdf),[],1);
        end
        [P T]=t2p(T);
    case {'Anders-Darl','AD'}  
        [YY I]=sort(Y);
        temp=((1:length(X)))./length(X);
        for i=size(Y,2):-1:1
            t_ecdf(I(:,i),i)=temp;
        end
        for i=size(perm,1):-1:1
            p_ecdf=NaN.*t_ecdf;
            for c=cc
                chi=find(X(perm(i,:))==c);
                [YY I]=sort(Y(chi,:));
                temp=(1:length(chi))./length(chi);
                for ii=size(Y,2):-1:1
                    p_ecdf(chi(I(:,ii)),ii)=temp;
                end
            end
            w=[t_ecdf.*(1-t_ecdf)].^.5;
            T(i,:)=sum(abs(p_ecdf-t_ecdf).*w);
        end
        [P T]=t2p(T);
end