%% MAIN SCRIPT FOR SPEARMAN BASED CORRELATIONS 
close all;clear all;
restoredefaultpath;
%% Load paths and data
addpath(genpath('/IMAGEN/AIDFUNCTIONS/'));
% addpath(genpath ('/usr/local/micapollo01/MIC/DATA/SHARED/STAFF/hhoske1/Students/svloem0/Software/AIDFUNCTIONS'))
% path to shared file:
% /usr/local/micapollo01/MIC/DATA/SHARED/STAFF/hhoske1/Students/svloem0/Software/AIDFUNCTIONS
studypath = '/usr/local/micapollo01/MIC/DATA/SHARED/STAFF/hhoske1/Students/svloem0/Syndromes/GWAS2/';
savepath = '/usr/local/micapollo01/MIC/DATA/STUDENTS/svloem0/Spearman_correlations/';
addpath('/usr/local/micapollo01/IMAGEN_DATA/SHARED/sgoova5/STUDENTS/Silke/functions/');
addpath('/usr/local/micapollo01/MIC/DATA/STUDENTS/svloem0/');
load('/usr/local/micapollo01/IMAGEN_DATA/SHARED/sgoova5/STUDENTS/Silke/SYNGWAS72.mat');
%% Load Syndromic GWAS data
% test = clearAllExcept(who,{'SYN'});for i=1:1:length(test),eval(['clear ' test{i}]);end
nSyndromes = 72; % 72
% if ~exist('SYN','var'), SYN = LoadMVPerSyndromeGWASmichp03(nSyndromes);end
% if ~exist('SYN','var'), SYN = LoadMVPerSyndromeGWAS_HH(nSyndromes);end %aangepast 
%% ALL SYNDROMES AT ONCE :)
in = load('/usr/local/micapollo01/IMAGEN_DATA/SHARED/hmatth5/SYNDROMETRAITGWAS2021/GWAS/SYNTRAITSMV4GWAS_17-Jan-2022.mat');
test = load('/usr/local/micapollo01/IMAGEN_DATA/SHARED/sgoova5/STUDENTS/Silke/LDblocks_EURO.mat');
hm3noMHC = readtable('/usr/local/micapollo01/IMAGEN_DATA/SHARED/sgoova5/STUDENTS/Silke/w_hm3.noMHC.snplist','FileType','text');
% load('/usr/local/micapollo01/IMAGEN_DATA/SHARED/sgoova5/STUDENTS/Silke/SYNGWAS72.mat')
% in = load([studypath 'Phenotypes/SYNTRAITSMV4GWAS.mat']); % Je moet eens
% checken of deze 72 syndromen bevat

GWAS = SYN;



% 1) Prepare P-values
ALLP = squeeze(SYN.P(:,:,:,3));
ALLP = max(ALLP,[],3);
ALLP = double(ALLP)/10000; % rows: SNPs – columns: traits

% 2) Match SNPs with hapmap set of SNPs
[~,ind21] = vlookupFast(SYN.RS,hm3noMHC.SNP);
CHR = SYN.CHR(ind21);
POS = SYN.POS(ind21);
ALLPhm3 = ALLP(ind21,:);

% 3) Run
[GC,pGC,seGC] = geneticCorrelationLDblocksWithSE(CHR,POS,ALLPhm3,ALLPhm3,test.REF.LDblocks,'mean');

% 4) Check
Names = in.syndromeNames;

imagesc(GC); colormap(parula); colorbar;
set(gcf,'Position',[400,100,1250,1200]);  set(gca, 'FontSize', 10);
xticks(1:nSyndromes); yticks(1:nSyndromes);
xticklabels(Names); yticklabels(Names); 

imagesc(-log10(pGC)); colormap(parula); colorbar;
set(gcf,'Position',[400,100,1250,1200]);  set(gca, 'FontSize', 10);
xticks(1:nSyndromes); yticks(1:nSyndromes);
xticklabels(Names); yticklabels(Names); 

% 5) Save
save('/usr/local/micapollo01/MIC/DATA/STAFF/sgoova5/tmp/GCSyn.mat','GC','pGC','seGC','Names');




