%% Investigating LEFT - RIGHT asymmetry
close all;clear;
%%
DATA_DIR = '../SAMPLE_DATA/';
RESULTS_DIR = '../results/demo_asymmetry/';
% THREADS = 8;
% samplesIndices = 1:1000;
nPicks = 10;
nSamplesPerPick = [200];
% nSamplesPerPick = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600];
% reduce = 0.05;
% reduce = 0.01;
% nRep = 3;
nIter = 1000;
THREADS = 8;
samplesIndices = 'all';
reduce = 0.1;
nRep = 10;
% Define following when nRep>1 -> no use of AMMI
% nIter = 5000;
%nSamplesPerPick = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600];
%nPicks = 10;
performExperiments = 0;
restoredefaultpath;
loadWhileInLab = 0;
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
genopath = [DATA_DIR, 'IMAGEN/BRAIN/UKBIOBANK/GENOTYPES/'];
covGenoPhenoPath = [DATA_DIR, 'IMAGEN/BRAIN/UKBIOBANK/COVARIATES/'];

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
    regphenopath = [phenopath Regions{r} '/'];
    disp(['PROCESSING BRAIN REGION: ' Regions{r}]);
    if ~loadWhileInLab
        DATA{r} = load([regphenopath 'STAGEOODATA']);
        DATA{r}.Region = DATA{r}.Region.Region;
    else
        DATA{r} = load([regphenopath 'STAGEOODATA']);
    end
end

Template = clone(DATA{1}.Region.AvgShape);

LH = DATA{1}.Region.AlignedShapes;
RH = DATA{2}.Region.AlignedShapes;
RH(:,1,:,:) = -1*RH(:,1,:,:);
%%
assert(min(size(LH,3), sum(str2double(DATA{1}.Region.IID) == str2double(DATA{2}.Region.IID))) == size(LH,3))
phenoIID = DATA{1}.Region.IID(1:size(LH,3));
%%

Render = cell(1,2);
for r=1:nR
    Render{r} = load([regphenopath 'RENDERMATERIAL.mat']);
end
brainSurface = Render{1};
%%
% Mirror template on y axis to express the right hemisphere instead of the left
brainSurface.RefScan.Vertices(:,1) = -1 * brainSurface.RefScan.Vertices(:,1);
[landmarksIndices, reducedFaces, toUpsampleLandmarksIndices]  = getDownsampledLandmarksIndices(brainSurface.RefScan,reduce,true);

%%
reducedLH = LH(landmarksIndices,:, :);
reducedRH = RH(landmarksIndices,:, :);
%%

%%
c= 1;
checkTemplate =  shape3D;
checkTemplate.Vertices = reducedLH(:,:,1);
checkTemplate.Faces = reducedFaces;
checkTemplate.ColorMode = "Indexed";
checkTemplate.Material = 'Dull';
checkTemplate.ViewMode = 'solid';
checkTemplate.Visible = true;
checkTemplate.PatchHandle.FaceColor = 'flat';
viewer = checkTemplate.viewer;

%%
clear DATA;
%%
if ischar(samplesIndices) && strcmp(samplesIndices, 'all')
    samplesIndices = 1:size(reducedLH,3);
end
nSamples = length(samplesIndices);
%%

%%


reducedTemplate= shape3D;
reducedTemplate.Vertices = brainSurface.RefScan.Vertices(landmarksIndices, :) ;
reducedTemplate.Faces = reducedFaces;
%%
[AlignedShapes,AvgShape,CentroidSizes] = GeneralizedProcrustesAnalysis(cat(3, reducedLH, reducedRH), reducedTemplate,3,true,false,true,false);

clear LH  RH;
%%
pdim = size(AlignedShapes,3)/2;

Ns = floor(1/reduce);
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
% display3DLandmarksArrows(Template, AvgShape);
%%
ReducedShapes = AlignedShapes(:, :, [ samplesIndices, (samplesIndices +size(AlignedShapes,3)/2)] );

reducedTemplateAdjacency = Template.Adjacency;
Shapes = permute(ReducedShapes,[2 1 3]);
%%
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

%%
pheno.diff = X2-X1;
pheno.iid = phenoIID;
save('pheno.mat', 'pheno');

%% TWO WAY PROCRUSTES ANOVA ON RANDOM SUBSETS OF THE DATA
if nRep == 1
    out = computeAmmiModel(ReducedShapes);
    nSamplesPerPick = nSamples;
else
    [setOut, avgOut,stdOut] = AsymmetryAnalysisOnSubsets(X1,X2,nSamplesPerPick,nPicks, nIter,mult,1);
    out = avgOut;
    
    if performExperiments
        [retS, retL, retR, retP] = ProcrustesAnova2WayAsymmetryDebuggingExperiments(X1, X2, mult); %#ok<UNRCH>
    end
end
%%
% Upsampling
outu = upsampleAnovaStats(out, reducedTemplateAdjacency, landmarksIndices);

showstruct = outu;
showPerm=1;
data = ProcrustesAnova2WayAsymmetryOutputProcess(...
    brainSurface, showstruct, nSamplesPerPick , showPerm, [RESULTS_DIR 'data_' experimentName '.mat']);

%%
rawF = outu.Raw.F;
%%
f = visualizeBrainAsymmetryData(data,[RESULTS_DIR 'results_' experimentName]);

%%
if nRep == 1
    phenoT = table(phenoIID, out.Raw.F');
    writetable(phenoT,[RESULTS_DIR 'fluctuatingAMMI.txt'],'WriteVariableNames',false,'Delimiter',' ');
end

%%
system('git add *');
message = ['AutoUpdate ' datestr(datetime('now'))];
system(['git commit -m "' message '"']);
system(['git push origin']);


