%% MAIN SCRIPT FOR SPEARMAN BASED CORRELATIONS 
%% DEVELOP GWAS INTERPRETATION SCRIPTS
close all;clear all;
restoredefaultpath;
addpath(genpath('/IMAGEN/AIDFUNCTIONS/'));
%% dedicated
test = clearAllExcept(who,{'SYN'});for i=1:1:length(test),eval(['clear ' test{i}]);end
%% Loadinf FACE GWAS, in our case the MV
nSyndromes = 61;
if ~exist('SYN','var'), SYN = LoadMVPerSyndromeGWASmichp03(nSyndromes);end
%% SYNDROOM SELECTEREN
in = load('/IMAGEN/FACE/EUROGWAS/SYNDROME_PHENOTYPES/SyndromeTraitMatrices29-May-2020.mat');
syndromeID = 13;
Name = in.syndromeNames{syndromeID};
GWAS = cell(1,2);
GWAS{1} = SYN;
GWAS{1}.P = squeeze(SYN.P(:,syndromeID,:,:));
GWAS{1}.Name = Name;
syndromeID = 6;
Name = in.syndromeNames{syndromeID};
GWAS{2} = SYN;
GWAS{2}.P = squeeze(SYN.P(:,syndromeID,:,:));
GWAS{2}.Name = Name;
%% LOADING REFERENCE GENOME DATA
REF = load('/IMAGEN/GENOMEREFMATERIALS/LDblocks_EURO.mat');% loads the definition of LD blocks
hm3noMHC = readtable('/IMAGEN/GENOMEREFMATERIALS/w_hm3.noMHC.snplist','FileType','text');% Determine overlap with BRAIN GWAS
factor = 10000; % needed to convert uint32 -log10(pvalues) to double
%% SELECTING THE BEST META P
factor = 10000;
for i=1:2
   GWAS{i}.P = squeeze(GWAS{i}.P(:,:,3));
   GWAS{i}.P = max(GWAS{i}.P,[],2);
   GWAS{i}.P = double(GWAS{i}.P)/factor;% convert to double from uint32 (not these are -log10(p-values) and not original p-values)
end
%% REDUCING GWAS to the list of SNPs from the Hapmap3 references without MHC clusters (also used as input to LDSC)
[~,ind21] = vlookupFast(GWAS{1}.RS,hm3noMHC.SNP);
for i=1:2
    GWAS{i}.CHR = GWAS{i}.CHR(ind21);
    GWAS{i}.POS = GWAS{i}.POS(ind21);
    GWAS{i}.RS = GWAS{i}.RS(ind21);
    GWAS{i}.A1 = GWAS{i}.A1(ind21);
    GWAS{i}.A2 = GWAS{i}.A2(ind21);
    GWAS{i}.MAF = GWAS{i}.MAF(ind21);
    GWAS{i}.N = GWAS{i}.N(ind21,:);
    GWAS{i}.P = GWAS{i}.P(ind21,:);
end
%% GETTING SPEARMAN CORRELATION
CHR = GWAS{1}.CHR;
POS = GWAS{1}.POS;
V1 = GWAS{1}.P;
V2 = GWAS{2}.P;
type = 'mean';
% NOTE V1 and V2 can be a block of p-values, so one collumn per syndrome,
% and if V1 = V2, you will run all pairwise correlations, look at
% implementation!
[GC,pGC,seGC] = geneticCorrelationLDblocksWithSE(CHR,POS,V1,V2,REF.LDblocks,type);
%% THE END



