function [P, T,options] = NP_ReM(Y,X,DES,B,tail,options)
%NONPARAMETRIC permutation test for repeated measures
%P = NP_ReM(Y,X,DES) performs m nonparametric permutation tests of the 
%null hypothesis of equality of the mean of the m response variables
%observed on a number of occasions, usually according to time, so that
%successive responses are dependent. Y is the [n*C]-by-m data matrix
%(n sample size, C number of repeated measures, m variables number) of
%responses. X is a vector of length n*C. The modality of X indicates
%to which time the corresponding row of Y belongs. If X is a scalar, then
%X indicates the number of repeated measures C. DES is the C-by-K 
%Design Matrix of comparisons (K is the number of comparisons). By default
%DES is the C-by-C-1 matrix:
%
%                            [ 1  0  0 ...
%                             -1  1  0 ...
%                              0 -1  1 ...]
%
%where 1 and -1 in any column indicate which measurements to compare. With
%the default setting, each variable at time t is compared with the same 
%variable at time t+1. If DES is a string, the following values are allowed:
%
%                  All: All comparisons
%                  Bal: Base line vs others (Base line is the first one)
%                  Seq: Sequential comparisons (1vs2, 2vs3, ...C-1vsC)
%                  Trd: Trend (1vs2...C, 12vs3...C, 123vs4...C)
%
%by default DES = 'Seq'.
%P is the (B+1)-by-m-by-K P-values matrix, where B is the number of random
%permutations set equal to 1000 by default. The last row of P is
%related to the observed sample. The default alternative of the test is 
%"means are not equal". We can change the default alternative by
%parameter TAIL (see below)
%
%P = NP_ReM(Y,X,DES,B) performs B random permutations.If B is a scalar,
%it indicates the number of random permutations. If B is a (B+1)-by-n matrix,
%it indicates the permutation sample space; in this case, each row is a
%(random) vector of 1 and -1. The last row must be the vector ones(1,n).
%if B is equal to 0 all possible permutations are done.
%
%P = NP_ReM(Y,X,DES,B,TAIL) the default alternative is "mean is not MU"
%With TAIL it is possible to specify other alternatives:
%
%               TAIL =  0,  alternative: "mean is not MU".
%               TAIL =  1,  alternative: "mean is greater than MU"
%               TAIL = -1,  alternative: "mean is less than MU"
%   
%where MU is the mean of the distribution (see options.MU below).
%If TAIL is an m-by-K matrix, the element(m,k) corresponds to variable m
%in the k comparison. If TAIL is an m length vector, each element
%corresponds to a variable (same tail for each comparison). If it is a
%scalar, this is considered for each variable.
%
%P = NP_ReM(...,'OPTIONS') Possible settings for options are:
%
%	options.OUT =  1 print the observed p-values, 0 don't show it. OUT = 1 by
%                  default.
%   options.Pobs = 1 compute p-values only for observed (not permuted)
%                  data. This save a considerable amount of time if the number
%                  of variables or permutations is high. Pobs = 0 by default
%                  (compute p-values for both observed and permuted data).
%   options.MU   = mean considered for the testing problem. By default, the
%                  null hypothesis tests the mean of distribution equal to 0.
%                  With options.MU it is possible to specify different values
%                  for the mean. If MU is a 1-by-m vector, each element 
%                  corresponds to a variable; if it is a scalar, this is
%                  considered for each variable.
%   options.labels.dims{2} = {'Var1', ......,'Varm'} to customize the
%                  labels of p-values in the output for the m variables Y.
%                  By default the labels are Y1, ..., Ym.
%   options.labels.dims{3} = {'Label first comparison', ...,'Label K-th
%                  comparison'} to customize the labels of p-values in the
%                  output for the K comparisons.
%
%[P,T] = NP_ReM(...) also returns the values of the test statistics
%in the (B+1)-by-m-by-K T matrix. The last row is related to the observed
%sample.
%
%[...,OPTIONS] =  NP_ReM(...) saves the options used for the analysis. 
%OPTIONS is a structural array with the following structure:
%
%       options.labels.dimslabel labels for the p-values of the matrix P:
%		         First dimension label 'Random Permutation'
%		         Second dimension label 'VARIABLES Y'
%		         Third dimension label 'COMPARISONS'
%	    options.labels.dims labels for any variable of any dimension of the
%                matrix P
%       options.OUT see above
%       options.Pobs see above
%	    options.tail see above
%	    options.p.raw    P-values
%
%The permutation methodology can handle missing data (whether missing at 
%random or not). Please indicate this as NaN.
%       
%Example: 
%We consider the data set GUARDA in the matlab file guardaex.mat.
%Below is the data set with 4 repeated measure in 7 subjects with 4
%variables:
%
%               paz	     Time  LTA   RTA    LTP   RTP
%               Paz. 11	 0	   3.3	 2.9    10.1  5.1      
%               Paz. 12	 0	   3.4	 3.7    4.4	  7.8
%               Paz. 13	 0	   4.7	 5	    8.8	  6.3
%               Paz. 14	 0	   5.3	 4.6	19.7  16.6
%               Paz. 15	 0	   5.9	 1.5	4.1	  4
%               Paz. 16	 0	   5.7	 3.4	4.9	  2.8
%               Paz. 17	 0	   8.1	 4.8	8.4	  8.5
%               Paz. 11	 1	   2	 1.1	17.2  8.4
%               Paz. 12	 1	   3.7	 2.2	10.3  6.8
%               Paz. 13	 1	   2.6	 2.3	13.6  3.9
%               Paz. 14	 1	   4.6	 4.8	9.4	  11.2
%               Paz. 15	 1	   4.4	 2.7	4.7	  4.3
%               Paz. 16	 1	   5.1	 2.6	8.7	  2
%               Paz. 17	 1	   3.2	 4.7	9.7	  6.4
%               Paz. 11	 2	   3.9	 4.3	11.6  6.5
%               Paz. 12	 2	   2.5	 1.9	7.3	  10.7
%               Paz. 13	 2	   1.9	 2.6	5.6	  3.8
%               Paz. 14	 2	   2.3	 2.7	24.8  24.6
%               Paz. 15	 2	   4.3	 1.8	6.9	  6.9
%               Paz. 16	 2	   1.4	 1.2	7.3	  6.5
%               Paz. 17	 2	   4.6	 4.4	10.7  6.5
%               Paz. 11	 3	   4.8	 1.7	7.2	  3.9
%               Paz. 12	 3	   3.7	 6.9	5.6	  4.4
%               Paz. 13	 3	   4.7	 3.7	3.3	  2.1
%               Paz. 14	 3	   4.9	 3.5	9.6	  8.4
%               Paz. 15	 3	   5.2	 2	    6.6	  5.7
%               Paz. 16	 3	   4.9	 4.9	4.6	  2.4
%               Paz. 17	 3	   3.4	 4.4	11.1  9.1
%
%
%To test the equality of means (means equal to 0) of all possible pairs of
%measurements and for each variable:
%
%load guardaex
%x = guardaex(:,2);
%y = guardaex(:,3:end);
%[P,T,options] = NP_ReM(y,x,'All');
%
%which produce the output:
%
%                       VARIABLES
%       COMPARISONS		Y1        	Y2        	Y3        	RTP        
%           1vs2		0.026973   	0.16284    	0.4016     	0.33866    
%           2vs3		0.41359    	0.76124    	0.98501    	0.091908   
%           3vs4		0.071928   	0.27073    	0.038961   	0.062937   
%           1vs3		0.023976   	0.10589    	0.088911   	0.18581    
%           2vs4		0.091908   	0.2997     	0.11089    	0.32368    
%           1vs4		0.50649    	0.7982     	0.37962    	0.18082      
%
%To customize the labels of Y variables, we must define the variable
%labels in the options:
%
%options.labels.dims{2} = {'LTA','RTA','LTP','RTP'};
%
%and repeat the test with the options (if options is used, it must
%also specify the B and TAIL parameters):
%
%[P,T,options] = NP_2s(y,x,1000,0,options);
%
%which produce the output: 
%
%                       VARIABLES
%       COMPARISONS		LTA        	RTA        	LTP        	RTP        
%           1vs2		0.026973   	0.16284    	0.4016     	0.33866    
%           2vs3		0.41359    	0.76124    	0.98501    	0.091908   
%           3vs4		0.071928   	0.27073    	0.038961   	0.062937   
%           1vs3		0.023976   	0.10589    	0.088911   	0.18581    
%           2vs4		0.091908   	0.2997     	0.11089    	0.32368    
%           1vs4		0.50649    	0.7982     	0.37962    	0.18082      
%
%These data are also reported in the Excel file guardaex.xls, with column
%labels showing the names of the variables. To perform the same analysis 
%with this file, we use the structural arrays D generated with xlsimport:  
%
%[D]=xlsimport('guardaex');
%[P,T,options] = NP_ReM(D(:,3:end),D(:,2),'All');
%
%This produce the output:
%
%                       VARIABLES
%       COMPARISONS		LTA        	RTA        	LTP        	RTP        
%           1vs2		0.026973   	0.16284    	0.4016     	0.33866    
%           2vs3		0.41359    	0.76124    	0.98501    	0.091908   
%           3vs4		0.071928   	0.27073    	0.038961   	0.062937   
%           1vs3		0.023976   	0.10589    	0.088911   	0.18581    
%           2vs4		0.091908   	0.2997     	0.11089    	0.32368    
%           1vs4		0.50649    	0.7982     	0.37962    	0.18082      
%
%It is possible to perform the same analysis, with the same output as above,
%using the variable labels of the columns of the Excel files as input.
%Beforehand it is necessary before to declare global the variable D:
%
%[D]=xlsimport('guardaex');
%reminD(D)
%[P,T,options] = NP_ReM({'LTA','RTA','LTP','RTP'},'Time','All');
%
%See also NP_2s, NP_Cs, NP_StOrd, NPC, XLSIMPORT.
%
%This function is part of NPClib.
%NPC Lib Matlab library (companion software of F. Pesarin, L. Salmaso: Permutation Tests for Complex Problems: theory, applications and software, John Wiley and sons).
%developed by Livio Finos; consulting team: Rosa Arboretti, Francesco Bertoluzzo, Stefano Bonnini, Chiara Brombin, Livio Corain, Luigi Salmaso and Aldo Solari.
%E-mail for technical questions: livio@stat.unipd.it (Livio Finos).


if nargin < 2
    error('Not enough iNPut argument')
elseif nargin == 2
   DES='Seq';
   B=1000;
end



if nargin==6
    [Y, options]=getDopts(Y,2,options);
else
    [Y, options]=getDopts(Y,2);
end
[X]=getDopts(X,3);

if isempty(DES),      options.DES=DES_ReM_std(length(unique(X)),'Seq'); 
elseif ischar(DES),       options.DES=DES_ReM_std(length(unique(X)),DES);
else options.DES=DES;
end


    clear DES
options=get_options(ones(1,size(Y,2),size(options.DES,1),size(X,3)),'ReM',options);


if isscalar(X)
    C=X;
    X=1:C;
    X=reshape(repmat(X,size(Y,1)./C,1),size(Y,1),1);
    id=1:size(Y,1);
else
    [X id]=sort(X);
end



Y=Y(id,:);

[no]=tabulate(X);
n=no(:,2);
no=no(:,1);
n=min(n);
C=length(no);
[m]=size(Y,2);
if not(length(X)./n==C)
    fprintf('\n sample size have to be equals!!')
    return;
end
Y=reshape(Y,[n C m]);


[C k]=size(options.DES);
if nargin <= 3
   tail=zeros(m,k);
   B=1000;
elseif nargin == 4
   tail=zeros(m,k);
end


temp=zeros(n,m,k);
for i=1:m
   temp(:,i,:)=Y(:,:,i)*options.DES;
end
temp(find(isnan(temp)))=0;
Y=temp;

clear temp
if length(options.MU)==1
    Y=Y-repmat(options.MU,[n m k]);
else
    Y=Y-repmat(reshape(options.MU,[1 m 1]),[n 1 k]);
end

if length(B)==1
    B=Space_sign(n,B);
end
B=B./n;
T=B*Y(:,1:m*k);
T=reshape(T,size(B,1),m,k);
options.B=size(B,1);

if prod(size(tail))==1
    tail=repmat(tail,[m k]);
elseif size(tail,1)==1
    tail=repmat(tail,[m 1]);
elseif size(tail,2)==1
    tail=repmat(tail,[1 k]);
end

[P T]=t2p(T,options.Pobs,tail);
options.p.raw=P(end,:,:);
NP_out('Testing repeated measures',options);