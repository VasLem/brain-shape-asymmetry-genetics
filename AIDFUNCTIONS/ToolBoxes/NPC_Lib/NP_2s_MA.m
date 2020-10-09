function [P, T, options] = NP_2s_MA(Y,X,B,stats,comb_funct,tail,options)
%Nonparametric permutation one-way MULTI-ASPECT ANOVA test for equality in
%distribution of 2 samples.
%P = NP_2s_MA(Y,X) performs a one-way ANOVA to compare the means of two 
%samples. Y is the n-by-m data matrix (n=n1+n2 sample size, 
%m variables number) of responses. X is a vector of length n. The two 
%modalities of X indicate to which sample the corrisponding row of Y belongs.
%If X is an n-by-q matrix, the same analysis is replicated for each column.
%If X or Y are Struct arrays, the data have a structure as outputted by
%xlsimport.m. 
%If X or Y are cell arrays, the cells are labels of the variables of X 
%and Y and variables must be declared global variables using reminD.m
%(and have a structure as in xlsimport.m). 
%P is the (B+1)-by-m-by-q P-values matrix, where B is the number of
%random permutations set equal to 1000 by default. The last row of P is
%related to the observed sample. The default alternative of the test is 
%"means are not equal". The default alternative can be changed by
%parameter TAIL (see below)
%
%P = NP_2S_MA(Y,X,B) performs B random permutations. If B is a scalar,
%it indicates the number of random permutations. If B is a (B+1)-by-n matrix,
%it indicates the permutation sample space; in this case, each row is a
%(random) permutation of the vector 1:n. The last row has to be the vector
%1:n itself if we wish  to make inference on the observed data. B = 1000
%by default. 
%
%P = NP_2s(Y,X,B,stats) to perform a MULTI-ASPECT one-way ANOVA test. stats
%is a 1-by-s cell array, where s is the number of considered aspects. You
%Possible settings for stats are:
%
%       't'                         for standard t-statistic
%       'Mann-Whitney',  (or 'MW')  for rank-based statistic
%       'Mean_grps'      (or 'M' )  for means comparisons
%       'Median_grps'    (or 'Me')  for medians comparisons
%       'Kolmog-Smirn'   (or 'KS')  for Kolmogorov-Smirnov statistic
%       'Anders-Darl'    (or 'AD'}  for Anderson-Darling statistic.
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
%P = NP_2s_MA(Y,X,B,stats,comb_funct) to specify the combining
%function used to combine the test statistics or p-values in a global
%p-value. See NPC for details.
%
%P = NP_2s_MA(Y,X,B,stats,comb_funct,tail). The default alternative is
%"means are not equal". With TAIL it is possible to specify other
%alternatives: 
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
%   options.Pobs = 1 compute p-values only for observed (not permuted)
%                  data. This save a considerable amount of time if the number
%                  of variables or permutations is high. Pobs = 0 by default
%                  (compute p-values for both observed and permuted data).
%   options.labels.dims{2} = {'Var1', ......,'Varm'} to customize the
%                  labels of p-values in the output for the m variables Y.
%                  By defoalt the labels are Y1, ..., Ym
%   options.labels.dims{3} = {'Var1', ......,'Varq'} to customize the
%                  labels of p-values in the output for the q variables X.
%                  By default the labels are X1, ..., Xq
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
%The permutation methodology can handle with missing data 
%(whether missing at random or not). Please indicate this as NaN.
%
%This function is part of NPClib.
%NPC Lib Matlab library (companion software of F. Pesarin, L. Salmaso: Permutation Tests for Complex Problems: theory, applications and software, John Wiley and sons).
%developed by Livio Finos; consulting team: Rosa Arboretti, Francesco Bertoluzzo, Stefano Bonnini, Chiara Brombin, Livio Corain, Luigi Salmaso and Aldo Solari.
%E-mail for technical questions: livio@stat.unipd.it (Livio Finos).

if nargin==7
    [Y, options]=getDopts(Y,2,options);
else
    [Y, options]=getDopts(Y,2);
end

[X, options]=getDopts(X,3,options);

    sizY=size(Y);
    sizX=size(X);
    options=get_options(ones([1 prod(sizY(2:end)) prod(sizX(2:end))]),'partial',options);

    
if nargin <= 1
    error('Not enough iNPut argument')
elseif nargin == 2
   B=1000;
   tail=zeros([prod(sizY(2:end)) prod(sizX(2:end))]);
elseif nargin == 3
   tail=zeros([1 prod(sizY(2:end)) prod(sizX(2:end))]);
end

if length(tail)==1
    sizY=size(Y);
    sizX=size(X);
    tail=tail.*ones([1 prod(sizY(2:end)) prod(sizX(2:end))]);
end

[origsizeY]=size(Y);
N=origsizeY(1);
m=prod(origsizeY(2:end));
Y=reshape(Y,N,m);  

miss=isfinite(Y);

if isscalar(B)
    [B] = Space_perm(N,B);
end
    mancanti=find(sum(miss)<size(Y,1));
    nonmanc=setdiff((1:size(Y,2)),mancanti);
    
    opts.OUT=0;
if not(isempty(mancanti))
    [Pm Tm]=NP_2s(Y(:,mancanti),X,B,tail(mancanti),opts);
end

Y=Y(:,nonmanc);
tail=tail(:,nonmanc);
[N m1]=size(Y);
P=zeros(size(B,1),m1,length(stats)+1);
T=P;
for j=1:length(stats)
    switch stats{j}(1)
        case {':'}
            stats{j}=stats{j}(2:end);
 %           ['[P(:,nonmanc,j) T(:,nonmanc,j)]=NP_ANOVA1(' stats{j} ',X,B,opts);']
            eval(['[P(:,:,j) T(:,:,j)]=NP_2s(' stats{j} ',X,B,tail,opts);'])
        otherwise
            [P(:,:,j) T(:,:,j)]=NP_2s_stat(Y,X,B,stats{j},tail);
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
function [P, T]=NP_2s_stat(Y,X,perm,stat,tail)
cc=unique(X)';
n=[(sum(X==cc(1))) (sum(X==cc(2)))];

switch stat
    case {'Median_grps','Me'}
            me=median(Y);
            for i=size(perm,1):-1:1
                T(i,:)=(median(Y(perm(i,X==cc(1)),:))-me);
            end
            T(:,tail(:)==1)=-T(:,tail(:)==1);
        T(:,tail(:)==0)=abs(T(:,tail(:)==0));          
        [P T]=t2p(T);
    case {'Mean_grps','M'}
            m=mean(Y);
            for i=size(perm,1):-1:1
                T(i,:)=(mean(Y(perm(i,X==cc(1)),:))-m);
            end        
            T(:,tail(:)==1)=-T(:,tail(:)==1);
            T(:,tail(:)==0)=abs(T(:,tail(:)==0));   
            [P T]=t2p(T);

    case {'Mann-Whitney','MW','t'}
        if (strcmp(stat,'Mann-Whitney')| strcmp(stat,'MW')),Y=R(Y);end
        displayopt='off';
        for i=size(perm,1):-1:1
            T(i,:)=ttest2m(Y,X(perm(i,:)),tail);
        end
        T(:,tail(:)==-1)=-T(:,tail(:)==-1);
        T(:,tail(:)==0)=abs(T(:,tail(:)==0));
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
            %tail
            for ii=size(Y,2):-1:1
            if not(tail(ii)==0)
                T(i,ii)=max(tail(ii).*(p_ecdf(:,ii)-t_ecdf(:,ii)));
            else
                T(i,ii)=max(abs(p_ecdf(:,ii)-t_ecdf(:,ii)));
            end
            end
                
        end
        [P T]=t2p(T);
    case {'Anders-Darl','AD'}  
        [YY I]=sort(Y);
        temp=((1:length(X))-.5)./length(X);
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

%%%%%%%%%%%%%%%%
function [t]=ttest2m(Y,X,tail)
t=zeros(1,size(Y,2));
gr=unique(X);

dfx = sum(X==gr(1)) - 1;
dfy = sum(X==gr(2)) - 1;
dfe  = dfx + dfy;
msx = dfx .* var(Y(X==gr(1),:));
msy = dfy .* var(Y(X==gr(2),:));

difference = mean(Y(X==gr(2),:)) - mean(Y(X==gr(1),:));
pooleds    = sqrt((msx + msy) .* (1/(dfx + 1) + 1/(dfy + 1)) ./ dfe);

t = difference ./ pooleds;
t(tail== 0)= abs(t(tail== 0)); 