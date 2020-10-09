%% DEVELOP GWAS INTERPRETATION SCRIPTS
close all;clear all;
restoredefaultpath;
addpath(genpath('/IMAGEN/AIDFUNCTIONS/'));
savepath = '/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/2020/SYNDROMEGWAS/DATA/TMPRESULTS/';
%% dedicated
test = clearAllExcept(who,{'REFSNP' 'SYN'});for i=1:1:length(test),eval(['clear ' test{i}]);end
%% Loading BRAIN GWAS SNP DATA, as referece SNP list towards LD computations
if ~exist('REFSNP','var'), REFSNP = LoadBRAINGWASmichp03('info');end
in = load('/IMAGEN/BRAIN/UKBIOBANK/PHENOTYPES/SWH/AUX_SWH_PA80_20200513.mat');
IID = in.DB.IID;
%% Loadinf FACE GWAS, in our case the MV
nSyndromes = 61;
if ~exist('SYN','var'), SYN = LoadMVPerSyndromeGWASmichp03(nSyndromes);end
%% SYNDROOM SELECTEREN
in = load('/IMAGEN/FACE/EUROGWAS/SYNDROME_PHENOTYPES/SyndromeTraitMatrices29-May-2020.mat');
syndromeID = 6;
Name = in.syndromeNames{syndromeID};
GWAS = SYN;
GWAS.P = squeeze(GWAS.P(:,syndromeID,:,:));
%% OVERLAP MET REF GENOME 
[ind12,ind21] = vlookupFast([uint32(GWAS.CHR) GWAS.POS],[uint32(REFSNP.CHR) REFSNP.POS],true);
GWAS.CHR = GWAS.CHR(ind21);
GWAS.POS = GWAS.POS(ind21);
GWAS.RS = GWAS.RS(ind21);
GWAS.A1 = GWAS.A1(ind21);
GWAS.A2 = GWAS.A2(ind21);
%GWAS.MAF = GWAS.MAF(ind21);
GWAS.N = GWAS.N(ind21);
GWAS.P = GWAS.P(ind21,:,:);
GWAS.JOBID = REFSNP.JOBID(ind12);
%% MIAMI PLOT GENEREREN
ftmp = figure;
cmap = colormap('redbluecmap');
mymap = [cmap(2,:);cmap(end-1,:)];
for i=1:12
    mymap = [mymap;cmap(2,:);cmap(end-1,:)];
end
close(ftmp);
% SIGNIFICANCE THRESHOLDS
pGW = 5e-8;
pSWB = pGW/40;
pSWF = pGW/40;
mode = 'other';
factor = 10000;

PUS = squeeze(GWAS.P(:,1,3));
PUK = squeeze(GWAS.P(:,2,3));
PUS = single(PUS)./factor;
PUK = single(PUK)./factor;

[f] = SYNMiamiplotUnion(GWAS.POS,GWAS.CHR,PUS,GWAS.POS,GWAS.CHR,PUK,pGW,pSWB,pSWF,mode,mymap(1:23,:),[]);
title(gca,Name);
print(f,'-dpng','-r300',[savepath '/MIAMIPLOT_' Name]);close all;
%% PEAK DETECTIE
%% Establishe thresholds
%% SETTING SIGNIFICANCE parameters
pGW = 5e-8;
pCrit = pGW;
pSUG = 5e-7;
pSW = pGW/40;
%% PRESELECTION OF SUGGESTIVE SIGNALS ONLY
P = squeeze(GWAS.P(:,:,3));% select US and UK Meta outcomes
P = single(P)/factor;
bestP = max(P,[],2);% keep the best from US or UK
GWAS.bestP = bestP;
index = find(bestP>=-log10(pSUG));%4619
GWASSUG = myReduceGWAS(GWAS,index);
GWASSUG.SNP = getBRAINGWASGENOTYPES(GWASSUG,IID);
save([savepath 'GWASSUG_' Name  '.mat'],'GWASSUG','-v7.3');
%% INITIAL DETECTION and PEAK CLUMPING
distmerge = 10e3;% anything within this distance will be clumped
distnet = 1e6;% anything within this distance will be checked for LD
cT = 1;% no phenotypic correlation taken into account
ldT = 0.01;% anythin in LD higher than this will be clumped
PEAKS = ExtractPeakByLocationAndLDAndPCorr(GWASSUG,pCrit,distnet,distmerge,ldT,cT);
disp(PEAKS.nPeak);% 37
BKPEAKS = PEAKS;
%% MERGING PEAKS
PEAKS = BKPEAKS;
distnet = 10e6;
PEAKS = MergePeakByLocationAndLD(PEAKS,distnet,ldT);
disp(PEAKS.nPeak);% 36
%% EXTRA REMOVING FOR ROBUSTNESS
%MAKE SURE TE REMOVE only below SW
minsupport = 2;
PEAKS = RemovePeakByMinSupportAndSign(PEAKS,pSW,minsupport);
disp(PEAKS.nPeak);% 35
disp(length(find(PEAKS.bestP>=-log10(pSW))));% 18
%% 
MainTab = table;
MainTab.RS = PEAKS.RS(:);
MainTab.CHR = PEAKS.CHR(:);
MainTab.POS = PEAKS.POS(:);
MainTab.A1 = PEAKS.A1(:);
MainTab.A2 = PEAKS.A2(:);
%MainTab.MAF = double(PEAKS.MAF(:))/100;
MainTab.N = PEAKS.N(:);
MainTab.SUPPORT = PEAKS.nSupport(:);
MainTab.LOG10P = PEAKS.bestP(:);
save([savepath 'PEAKS_' Name '.mat'],'PEAKS','MainTab');
writetable(MainTab,[savepath 'PEAK_' Name '.xlsx'],'FileType','spreadsheet','WriteVariableNames',true,'Sheet',1);
%%









