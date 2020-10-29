function [adj_P, p_glob, options] = NPC_FWE(T,comb_funct,options)
%Nonparametric Combination Methodology with control of multiplicity.
%adj_P = NPC_FWE(T) performs a closed testing procedure controlling the
%Familywise Error Rate (FWE, the probability of rejecting any true null
%hypothesis contained in a subset of true null hypotheses). 
%By default it makes use of Tippet's combining function. adj_P is a 1-by-m
%vector of adjusted p-values where m is the number of partial hypotheses
%to be tested. T is the (B+1)-by-m p-values or statistics matrix. The last
%row is related to the observed sample. P must be a p-values or statistics
%matrix in accordance with the kind of combining function we use (see
%below)
%
%adj_P = NPC_FWE(P,comb_funct) comb_funct specify the combination function
%to use. Possible settings for comb_funct are:
%
%           'Direct'        (or 'D')  combining function of t-statistic
%                                     (Closed testing),
%           'Max-t'         (or 'M')  combining function of t-statistic
%                                     (step-down procedure),
%           'Fisher'        (or 'F')  combining function of p-values
%                                     (Closed testing),
%           'Liptak'        (or 'L')  combining function of p-values
%                                     (Closed testing),
%           'LiptakLog'     (or 'LL') combining function of p-values
%                                     (Closed testing),
%           'Tippett'       (or 'T')  combining function of p-values 
%                                     (min-p step-down procedure).
%           'Bonf-Holm'     (or 'B')  Bonferroni-Holm step-down procedure
%           'Fisher asympt' (or 'Fa') asymptotic Fisher combination:
%                                     -2sum(log(p_i);1=1,...k)~Chi2(2k).
%
%By default Tippet's combining function is used.
%
%%adj_P = NPC_FWE(...,'OPTIONS') Possible settings for options are:
%
%   options.OUT =  1 print the observed p-values, 0 don't show it. OUT = 1
%                  by default.
%   options.w =    1-by-m vector of weights for the combining function.
%                  w = [1 1 ...].*1/m by default. 
%   options.threshold = is a scalar ]0;1]. The combining function does not 
%                       consider values greater than OPTIONS.threshold 
%                       (it set the p-value at 1). For 'Direct' or 'Max-t'
%                       combining functions it is a scalar ]0;inf[. The
%                       combining function does not consider values below
%                       OPTIONS.threshold (it sets the p-value at 0). 
%                       For the 'Tippett' (/'Max-t') combining function the
%                       procedure stops when the first supra(sub)
%                       threshold is reached. This can lead to a
%                       considerable reduction in computation time. This 
%                       options does not affect 'Bonf-Holm' and 'Fisher
%                       asympt' procedures.
%       options.Combdims dimension of the combination. By default it is set
%                equal to ndims(T), i.e. the last dimension of 
%                matrix P.
%
%[adj_P p_glob] = NPC_FWE(...) also returns the global p-value.
%
%[...,OPTIONS] = NPC_FWE(...) saves the options used for the analysis.
%OPTIONS is a structural array with the following structure:
%
%       options.labels.dimslabel label for the dimensions of the P-value
%                 matrix P:
%                 First dimension label 'Random Permutation'
%                 Second dimension label 'VARIABLES Y'
%                 Third dimension label 'VARIABLES X'
%                 Fourth dimension label 'STRATA'
%       options.labels.dims labels for any variable of any dimension of the
%                P-value matrix P
%       options.OUT see above
%       options.Pobs see above
%       options.Combdims dimension of the combination. By default it is set
%                equal to ndims(T), i.e. the last dimension of 
%                matrix P.
%       options.p.raw    P-values
%
%Example:
%We consider the data set westfall_wallfingerNPCFWE in the matlab file 
%westfall_wallfingerNPCFWE.mat. Below is the data set with two samples (X = 1
%and X = 2) and three variables (Y1, Y2, and Y3):
%
%                   X       Y1      Y2      Y3
%                   0       14,4	7       4,3
%                   0       14,6	7,09	3,88
%                   0       13,8	7,06	5,34
%                   0       10,1	4,26	4,26
%                   0       11,1	5,49	4,52
%                   0       12,4	6,13	5,69
%                   0       12,7	6,69	4,45
%                   1       11,8	5,44	3,94
%                   1       18,3	1,28	0,67
%                   1       18      1,5     0,67
%                   1       20,8	1,51	0,72
%                   1       18,3	1,14	0,67
%                   1       14,8	2,74	0,67
%                   1       13,8	7,08	3,43
%                   1       11,5	6,37	5,64
%                   1       10,9	6,26	3,47
%
%First we perform a two-samples test with the NP_2s function:
%
%load westfall_wolfingerNPCFWE;
%y = westfall_wolfingerNPCFWE(:,2:4);
%x = westfall_wolfingerNPCFWE(:,1);
%P = NP_2s(y,x,1000,0);
%
%which produce the output:
%  
%________________________________________________________________
% Testing equality in distribution for 2 independent samples 
% for continuous (or dichotomous) variables (with or without missing values) 
% Number of Conditional Montecarlo: 1001
%            			            
%	VARIABLES Y 	Y1          	Y2          	Y3          
%	            	0.1139      	0.03497     	0.01399     
%
%Now we can combine the three partial p-values corrected
%for multiplicity with NPC_FWE to obtain a global test:
%
%adj_P = NPC_FWE(P,'F');
%
%which produce the output:
%
%________________________________________________________________
%Non Parametric Closed Testing procedure for the strong control of the
%Familywise Error Rate
%Number of Conditional Montecarlo: 1001
%            			            
%	VARIABLES Y 	             
%	Y1          
%	p-value     	0.11389      
%	   p-FWE    	 0.11389     
%	Y2          
%	p-value     	0.034965     
%	   p-FWE    	 0.052947    
%	Y3          
%	p-value     	0.013986     
%	   p-FWE    	 0.038961    
%
% Comb Funct		 Fisher
% p-GLOBAL	
%	            	            
%	            	0.031968  
%
%This function is part of NPClib.
%NPC Lib Matlab library (companion software of F. Pesarin, L. Salmaso: Permutation Tests for Complex Problems: theory, applications and software, John Wiley and sons).
%developed by Livio Finos; consulting team: Rosa Arboretti, Francesco Bertoluzzo, Stefano Bonnini, Chiara Brombin, Livio Corain, Luigi Salmaso and Aldo Solari.
%E-mail for technical questions: livio@stat.unipd.it (Livio Finos).

if nargin==1
    comb_funct='T';
end
if nargin<3
    options=get_options(T,'NPC');
else
    options=get_options(T,'NPC',options);
end

dims=size(T);
temp=ones(1,options.Combdims);
temp(1:length(dims))=dims;
dims=temp;
if isfield(options,'Combdims')
    T=split2dims2(T,options.Combdims);
    T=reshape(T,[dims(1) dims(options.Combdims) prod(dims(setdiff(2:length(dims),options.Combdims)))]);
else
    T=reshape(T,[dims(1:2) prod(dims(3:length(dims)))]);
end

[B m ntimes]=size(T);

TT=T;
options.p.raw=ones(1,m,0);
options.p.adj=ones(1,m,0);
options.p.glb=ones(1,1,0);

for it=1:ntimes
    T=TT(:,:,it);
%assegna raw p-values, threshold
switch comb_funct
    case {'Direct','D','d','Max-t','M','m'}
        P.raw=t2p(T,1);
        if isfield(options,'threshold'),
            T(T<options.threshold)=0;
        end
    case {'Bonf-Holm','B','b'}
        P.raw=T(:)';
        if isfield(options,'w')==0,     
            options.w=ones(1,length(T)); 
        end
    case {'Fisher asympt','Fa','fa','FA'}
        P.raw=T(:)';
        if isfield(options,'w'),     
            fprintf('\nWARNING: weigths are not considered using Fisher asymptotic combination')
            options=rmfield(options,'w');
        end
    otherwise
        P.raw=T(end,:);
        if isfield(options,'threshold'),
            T(T>options.threshold)=1;
        end
end    


        
%trasformazione di p-value e "sistemazione" dei dati
switch comb_funct
    case {'Direct','D','d'}
        CF='Direct';
    case {'Max-t','M','m'}
        CF='Max-t';
        [null ind]=sort(-T(B,:));
        clear null;
    case {'Fisher','F','f'}
        CF='Fisher';
        T=(T.*B -.5)/B;
        T=-log(T);
    case {'Liptak','L','l'}
        CF='Liptak';
        T=(T.*B -.5)/B;
        T=-norminv(T);
    case {'LiptakLog','LL','ll'}
        CF='LiptakLog';
        T=(T.*B -.5)/B;
        T=log((1-T)./T);
    case {'Tippett','T','t'}
        CF='Tippett';
        T=-T;
        if isfield(options,'w'),
            options.w=options.w.^-1;
        end
        [null ind]=sort(-T(B,:));
        clear null;
    case {'Bonf-Holm','B','b'}
        CF='Bonferroni-Holm';
        if m==1
            m=B;
            T=T';
        else
            T=T(B,:);
        end
    case {'Fisher asympt', 'Fa','fa','FA'}
        CF='Fisher (asymptotic)';
        if m==1
            m=B;
            T=-2*log(T)';
        else
            T=-2*log(T(B,:));
        end
end


%assegna pesi
if isfield(options,'w'),
    if var(options.w)>0
        T=T*diag(options.w);
    end
end

%closed testing
switch comb_funct
    case {'Fisher','F','f','Liptak','L','l','Direct','D','d','LiptakLog','LL','ll'}
        
        try 
            in_tst=ff2n(m);
        p_cmb_f=transpose(in_tst*T');
        catch
        error('\nSorry, the maximum variable size allowed by NPC_FWE with Fisher, Liptak, Direct or Liptaklog is exceeded. \nYou are allowed to use Tippett, Max-T, Bonf-Holm or Fisher asympt combining function (see the help)\n','%s')
        break
        end
        cmb_pvals=t2p(p_cmb_f,1);
        clear p_cmb_f;
        for i=1:m
            adj_P(i)=max(cmb_pvals(find(in_tst(:,i))));
        end
        p_glob=cmb_pvals(end);
    case {'Tippett','T','t','Max-t','M','m'}
        
        p_cmb_hoch=ones(B,m);
        
        if isfield(options,'threshold')
            switch comb_funct
                case {'Max-t','M','m'}
                    halt=sum(P.raw>=options.threshold);    
                case {'Tippett','T','t'}
                    halt=sum(P.raw<=options.threshold);
            end
        else
            halt=length(P.raw);
        end
        
        for i=1:halt
           p_cmb_hoch(:,i)=max(T(:,ind(i:m)),[],2);
        end
        adj_P(1,ind(1:halt))=t2p(p_cmb_hoch(:,1:halt),1);
        adj_P(1,ind(halt+1:end))=1;
        for i=2:halt
            adj_P(1,ind(i))=max([adj_P(1,ind(i)) adj_P(1,ind(i-1))],[],2);
        end
        p_glob=adj_P(1,ind(1));
    case {'Bonf-Holm','B','b'}
         [null ind]=sort(T);
         options.w=options.w(ind)./(sum(options.w(:))).*m;
         I(m:-1:1) = (cumsum(options.w(end:-1:1)))';
         null = null./options.w.*I;
        adj_P(ind(1))=null(1);
        for i=2:length(ind)
            adj_P(ind(i))=max(adj_P(ind(i-1)),null(i));
        end
        adj_P=min(adj_P,ones(1,m));
        p_glob=adj_P(1,ind(1));
    case {'Fisher asympt','Fa','fa','FA'}
        if m==1
            m=B;
            T=T';
        else
            T=T(B,:);
        end
        in_tst=ff2n(m);
        p_cmb_f=1-chi2cdf(in_tst*T',2.*sum(in_tst,2));
        adj_P=max(in_tst.*repmat(p_cmb_f,1,m),[],1);
        p_glob=p_cmb_f(end);
end

options.p.raw(:,:,it)=P.raw(:);
options.p.adj(:,:,it)=adj_P;
options.p.glb(:,:,it)=p_glob;
end

if isfield(options,'Combdims')
        options.p.raw=split2dims2(options.p.raw,options.Combdims);
        options.p.raw=reshape(options.p.raw,[1 dims(2:end)]);
        options.p.adj=split2dims2(options.p.adj,options.Combdims);
        options.p.adj=reshape(options.p.adj,[1 dims(2:end)]);
        options.p.glb=split2dims2(options.p.glb,options.Combdims);
        options.p.glb=reshape(options.p.glb,[1 dims(2:options.Combdims-1) 1 dims(options.Combdims+1:end)]);
end
    options.B=B;


 options.title='Non Parametric Closed Testing procedure for the strong control of the Familywise Error Rate';
    options.CF=CF;

NPC_out(options);

adj_P=options.p.adj;
