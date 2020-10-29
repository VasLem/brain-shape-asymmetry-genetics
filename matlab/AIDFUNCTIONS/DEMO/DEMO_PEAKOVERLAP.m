%% DEVELOP GWAS INTERPRETATION SCRIPTS
%close all;clear all;
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
A = load('/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/2020/SYNDROMEGWAS/DATA/TMPRESULTS/PEAKS_Cleft Lip Palate.mat');
tmp = load('/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/2020/SYNDROMEGWAS/DATA/TMPRESULTS/GWASSUG_Cleft Lip Palate.mat');
A.GWASSUG = tmp.GWASSUG;

B = load('/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/2020/SYNDROMEGWAS/DATA/TMPRESULTS/PEAKS_Achondroplasia.mat');
tmp = load('/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/2020/SYNDROMEGWAS/DATA/TMPRESULTS/GWASSUG_Achondroplasia.mat');
B.GWASSUG = tmp.GWASSUG;

%% setting overlap parameters
distnet = 1e6;
ldT = .2;
mindist = 10e3;
pGW = 5e-8;
pSUG = 5e-7;
%% A->B
out = overlapPeaksInGWASWithLD(A.PEAKS,B.GWASSUG,distnet,mindist,ldT);
sum(out.exactOverlap>=-log10(pGW))
sum(out.exactOverlap>=-log10(pSUG))
sum(out.proxOverlap>=-log10(pGW))
sum(out.proxOverlap>=-log10(pSUG))
Tab = table;
Tab.ExactSNP_LOG10P = out.exactOverlap;
Tab.ProxySNP_LOG10P = out.proxOverlap;
rs = cell(PEAKS.nPeak,1);
pos = nan*zeros(PEAKS.nPeak,1);
r2 = nan*zeros(PEAKS.nPeak,1);
for i=1:PEAKS.nPeak
    if isempty(out.proxInfo{i}), continue; end
    rs{i} = out.proxInfo{i}.RS{1};
    pos(i) = out.proxInfo{i}.POS;
    r2(i) = out.proxInfo{i}.R2;
end
Tab.Proxy_RS = rs(:);
Tab.Proxy_POS = pos(:);
Tab.R2 = r2(:);
OVERLAP = out;
OVERLAPTAB = [A.MainTab Tab];
outAB = out;
%save([savepath  '_OVERLAPINFACEGWAS.mat'],'PEAKS','OVERLAP','OVERLAPTAB');
%writetable(OVERLAPTAB,[outputpath prefix '_OVERLAPINFACEGWAS.xlsx'],'FileType','spreadsheet','WriteVariableNames',true,'Sheet',1);
%% A->B
out = overlapPeaksInGWASWithLD(B.PEAKS,A.GWASSUG,distnet,mindist,ldT);
sum(out.exactOverlap>=-log10(pGW))
sum(out.exactOverlap>=-log10(pSUG))
sum(out.proxOverlap>=-log10(pGW))
sum(out.proxOverlap>=-log10(pSUG))
Tab = table;
Tab.ExactSNP_LOG10P = out.exactOverlap;
Tab.ProxySNP_LOG10P = out.proxOverlap;
rs = cell(PEAKS.nPeak,1);
pos = nan*zeros(PEAKS.nPeak,1);
r2 = nan*zeros(PEAKS.nPeak,1);
for i=1:PEAKS.nPeak
    if isempty(out.proxInfo{i}), continue; end
    rs{i} = out.proxInfo{i}.RS{1};
    pos(i) = out.proxInfo{i}.POS;
    r2(i) = out.proxInfo{i}.R2;
end
Tab.Proxy_RS = rs(:);
Tab.Proxy_POS = pos(:);
Tab.R2 = r2(:);
OVERLAP = out;
OVERLAPTAB = [B.MainTab Tab];
outBA = out;
%save([savepath  '_OVERLAPINFACEGWAS.mat'],'PEAKS','OVERLAP','OVERLAPTAB');
%writetable(OVERLAPTAB,[outputpath prefix '_OVERLAPINFACEGWAS.xlsx'],'FileType','spreadsheet','WriteVariableNames',true,'Sheet',1);
%%
out = outAB;
PEAKS = A.PEAKS;
index = find(out.proxOverlap>=-log10(pSUG));
PEAKS.RS = PEAKS.RS(index);
PEAKS.POS = PEAKS.POS(index);
PEAKS.CHR = PEAKS.CHR(index);
PEAKS.A1 = PEAKS.A1(index);
PEAKS.A2 = PEAKS.A2(index);
PEAKS.N = PEAKS.N(index);
PEAKS.bestP = PEAKS.bestP(index);
PEAKS.P = PEAKS.P(index,:,:);
PEAKS.SNP = PEAKS.SNP(index,:);
PEAKS.Support = PEAKS.Support(index);
PEAKS.nSupport = PEAKS.nSupport(index);
PEAKS.nPeak = length(index);
PEAKSA = PEAKS;
out = outBA;
PEAKS = B.PEAKS;
index = find(out.proxOverlap>=-log10(pSUG));
PEAKS.RS = PEAKS.RS(index);
PEAKS.POS = PEAKS.POS(index);
PEAKS.CHR = PEAKS.CHR(index);
PEAKS.A1 = PEAKS.A1(index);
PEAKS.A2 = PEAKS.A2(index);
PEAKS.N = PEAKS.N(index);
PEAKS.bestP = PEAKS.bestP(index);
PEAKS.P = PEAKS.P(index,:,:);
PEAKS.SNP = PEAKS.SNP(index,:);
PEAKS.Support = PEAKS.Support(index);
PEAKS.nSupport = PEAKS.nSupport(index);
PEAKS.nPeak = length(index);
%%
PEAKS.RS = [PEAKS.RS PEAKSA.RS];
PEAKS.POS = [PEAKS.POS PEAKSA.POS];
PEAKS.CHR = [PEAKS.CHR PEAKSA.CHR];
PEAKS.A1 = [PEAKS.A1 PEAKSA.A1];
PEAKS.A2 = [PEAKS.A2 PEAKSA.A2];
PEAKS.N = [PEAKS.N PEAKSA.N];
PEAKS.bestP = [PEAKS.bestP PEAKSA.bestP];
PEAKS.P = cat(1,PEAKS.P,PEAKSA.P);
PEAKS.SNP = [PEAKS.SNP; PEAKSA.SNP];
PEAKS.Support = [PEAKS.Support PEAKSA.Support];
PEAKS.nSupport = [PEAKS.nSupport PEAKSA.nSupport];
PEAKS.nPeak = length(PEAKS.RS);

%%
PEAKS = MergePeakByLocationAndLD(PEAKS,distnet,ldT);
peaks = PEAKS;% NOTE PEAKS ARE UNSORTED 
%% sorting peaks
% PEAKS = [];
% PEAKS.RS = [];
% PEAKS.CHR = [];
% PEAKS.POS = [];
% PEAKS.SNP = [];
% %PEAKS.source = [];
% for c=1:23
%     % c=1;
%     index = find(peaks.CHR==uint8(c));
%     if isempty(index), continue; end
%     [~,ind] = sort(peaks.POS(index),'ascend');
%     PEAKS.RS = [PEAKS.RS(:); peaks.RS(index(ind))];
%     PEAKS.CHR = [PEAKS.CHR(:); peaks.CHR(index(ind))];
%     PEAKS.POS = [PEAKS.POS; peaks.POS(index(ind))];
%     PEAKS.SNP = [PEAKS.SNP; peaks.SNP(index(ind),:)];
%     %PEAKS.source = [PEAKS.source;peaks.source(index(ind))];
% end







