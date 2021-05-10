%% Investigating LEFT - RIGHT asymmetry
close all;clear;
%%
DATA_DIR = '../SAMPLE_DATA/';
% THREADS = 8;
% samplesIndices = 1:1000;
% nPicks = 10;
% nSamplesPerPick = [200];
% nSamplesPerPick = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600];
% reduce = 0.05;
% reduce = 0.01;
% nRep = 3;
% nIter = 1000;
THREADS = 20;
samplesIndices = 1:100;
reduce = 0.1;
nRep = 1;
% Define following when nRep>1 -> no use of AMMI
% nIter = 5000;
%nSamplesPerPick = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600];
%nPicks = 10;
% performExperiments = 0;
restoredefaultpath;
nSamples = length(samplesIndices);
loadWhileInLab = 1;
addpath(genpath('AIDFUNCTIONS'));
% addpath(genpath('/IMAGEN/AIDFUNCTIONS/'));

%% SETTING UP COMPUTATION POWER
try
    parpool('local',THREADS);
    %parpool('local',8);
catch
end
%% GETTING SOME INFO ON THE BRAIN TEMPLATE


in = load([DATA_DIR, 'IMAGEN/BRAIN/HumanConnectomeProject/SubcorticalMask_HCP.mat']);

phenopath = [DATA_DIR, 'IMAGEN/BRAIN/UKBIOBANK/PHENOTYPES/'];

MASK = in.index;
nVertices = length(MASK);
Regions = {'LH' 'RH'};
nR = length(Regions);


%% Extract mat file of serialized objects
% strDATA = cell(1,2);
% strRender = cell(1, 2);
% savePath = [DATA_DIR 'asym_data.mat'];
% disp(['Extracting mat file of serialized objects to ' savePath]);
% for r=1:nR
%     %r=2
%     regphenopath = [phenopath Regions{r} '/'];
%     strDATA{r} =  load([regphenopath 'STAGE00DATA']);
%     strDATA{r}.Region.AvgShape = strDATA{r}.Region.AvgShape.obj2struc();
%     strRender{r} = load([regphenopath 'RENDERMATERIAL.mat']);
%     strRender{r}.RefScan = strRender{r}.RefScan.obj2struc();
%
% end
% save(savePath,  'strDATA', 'strRender');
% clear strDATA strRender;

%% WITH THE UNCORRECTED DATA SET

% phenopath = '/IMAGEN/BRAIN/UKBIOBANK/PHENOTYPES/';


DATA = cell(1,2);

for r=1:nR
    %r=2
    disp(['PROCESSING BRAIN REGION: ' Regions{r}]);
    regphenopath = [phenopath Regions{r} '/'];
    DATA{r} = load([regphenopath 'STAGE00DATA']);
end

Template = clone(DATA{1}.Region.AvgShape);
LH = DATA{1}.Region.AlignedShapes;
RH = DATA{2}.Region.AlignedShapes;

RH(:,1,:,:) = -1*RH(:,1,:,:);





Render = cell(1,2);
for r=1:nR
Render{r} = load([regphenopath 'RENDERMATERIAL.mat']);
end
brainSurface = Render{1};
%%





%% Subsampling space to match R version which has memory issues to pick the whole space
clear DATA;
if loadWhileInLab
 load([phenopath 'intAligned.mat'],'mult','intLHAligned','intRHAligned');
 LH = double(intLHAligned )* mult;
 RH = double(intRHAligned) * mult;
 clear intLHAligned intRHAligned;
 AlignedShapes = cat(3,LH,RH);
end
[AlignedShapes,AvgShape,CentroidSizes] = GeneralizedProcrustesAnalysis(cat(3, LH, RH), Template,3,true,false,true,false);

clear LH  RH;
%%
pdim = size(AlignedShapes,3)/2;

Ns = floor(1/reduce);
landmarksIndices = 1:Ns:pdim;
nLandmarks = length(landmarksIndices);

if nRep > 1
    experimentName = [num2str(length(samplesIndices)) '_' num2str(Ns) '_' ...
    num2str(nIter)  '_' num2str(nPicks) '_' num2str(length(nSamplesPerPick)) '_' num2str(nSamplesPerPick(1))...
    '_' num2str(nSamplesPerPick (end)) '_' num2str(nRep)];
else
    experimentName = [num2str(length(samplesIndices)) '_' num2str(Ns) ...
    '_' num2str(nRep)];
end
%%






%%
% display3DLandmarksArrows(Template, AvgShape);
%%
ReducedShapes = AlignedShapes(landmarksIndices, :, [ samplesIndices, (samplesIndices +size(AlignedShapes,3)/2)] );

reducedTemplateAdjacency = Template.Adjacency;
Shapes = permute(ReducedShapes,[2 1 3]);
Shapes = reshape(Shapes,size(Shapes,1)*size(Shapes,2),size(Shapes,3))';
%%
%%
% mag = var(Shapes,0,2);
if nRep > 1
    percent_difference_test_retest = load('../SAMPLE_DATA/MeasError/percent_difference_test_retest');
    percent_difference_test_retest = round(percent_difference_test_retest.percent_difference_test_retest,2);
    percent_difference_test_retest = percent_difference_test_retest(landmarksIndices);
    percent_difference_test_retest = repmat(percent_difference_test_retest,1,3)';
    percent_difference_test_retest = percent_difference_test_retest(:)';
    RepShapes = zeros(size(Shapes,1),size(Shapes,2),nRep,'single');
    for i=1:1:nRep
        RepShapes(:,:,i) = single(Shapes) +single(randn(size(Shapes,1),size(Shapes,2)).*single(Shapes).*percent_difference_test_retest/100);
%         RepShapes(:,:,i) = single(Shapes) +single(randn(size(Shapes,1),size(Shapes,2)).*0.2.*mag);
        
    end
else
    RepShapes = zeros(size(Shapes,1),size(Shapes,2),1,'single');
    RepShapes(:,:,1) = Shapes;
end
mult = double(intmax('int16')) / (max(abs(RepShapes),[],'all'));
RepShapesInt16 = int16(RepShapes.*mult);
if nRep > 1
    
figure; histogram(reshape(RepShapesInt16 - int16(mult * Shapes), 1,[])); title({'Landmark Coordinate dislocation','for generated replications'});
end

%%
clear Shapes;
clear RepShapes;
%%
X1 = RepShapesInt16(1:nSamples,:,:);
X2 = RepShapesInt16(nSamples+1:end,:,:);



%% TWO WAY PROCRUSTES ANOVA ON RANDOM SUBSETS OF THE DATA
if nRep == 1
   out = computeAmmiModel(ReducedShapes);

else
    [setOut, avgOut,stdOut] = AsymmetryAnalysisOnSubsets(X1,X2,nSamplesPerPick,nPicks, nIter,mult,1);
    out = avgOut;

    if performExperiments
    [retS, retL, retR, retP] = ProcrustesAnova2WayAsymmetryDebuggingExperiments(X1, X2, mult); %#ok<UNRCH>
    end
end
% Upsampling
outu = upsampleAnovaStats(out, reducedTemplateAdjacency, landmarksIndices);

showstruct = outu;
showPerm=1;
data = ProcrustesAnova2WayAsymmetryOutputProcess(...
    brainSurface, showstruct, nSamplesPerPick , showPerm, ['../results/demo_asymmetry/data_' experimentName '.mat']);

%%
rawF = out.Raw.F;
save('rawF.mat','rawF');
system('git add rawF.mat')
system('git add -u');
message = ['AutoUpdate ' datestr(datetime('now'))];
system(['git commit -m "' message '"']);
system(['git push origin']);
%%
f = visualizeBrainAsymmetryData(data,['../results/demo_asymmetry/results_' experimentName]);

%%

genotypes = readtable("/home/vaslem/code/imagen/SAMPLE_DATA/PLINKPRUNEDMERGED/ukb_sel19908.fam",'FileType','text');
%%i
if nRep==1
ids = genotypes.Var1(1:1000);
T = table(ids, out.Raw.F');
writetable(T,"fluctuatingAMMI.txt",'WriteVariableNames',false,'Delimiter',' ')
end
%%
system('git add *');
message = ['AutoUpdate ' datestr(datetime('now'))];
system(['git commit -m "' message '"']);
system(['git push origin']);


