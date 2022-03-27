%ENV TO SET:
%%%%%%%%%%%%
%DATA_ROOT: the directory of the DATA (../SAMPLE_DATA)
%DATASET_INDEX: dataset to use, 1 or 2 (1)
%RESULTS_ROOT: the directory of the results (../results)
%SCRATCH_ROOT: the directory to use for the intermediate files (../results)
%THREADS: Number of threads to use (auto set by local)
%%%%%%%%%%%%
%% Applying Hierarchical Clustering to Brain Symmetry Related information
close all;clear;
if ~isdeployed
    restoredefaultpath;
    addpath(genpath('AIDFUNCTIONS'));
    addpath(genpath('FUNCTIONS'));
    addpath(genpath('.'));
end
NUM_LEVELS = 4;
MAX_NUM_PCS = 500;
DEFAULT_DATASET_INDEX = 1;
REDUCTION_RATE = 1;
DEFAULT_COMPONENT = 'asymmetry'; %asymmetry,symmetry
SEED = 42;

rng(SEED); % For reproducible results
if endsWith(cd, "AIDFUNCTIONS/DEMO")
    cd("../..")
end

DATA_DIR = '../SAMPLE_DATA/';
DATASET_INDEX = DEFAULT_DATASET_INDEX;

switch DATASET_INDEX
    case 1
        UKBIOBANK = 'UKBIOBANK';
        DATASET_NAME = 'STAGE00DATA';
        GENO_ID = 'sel19908';
    case 2
        UKBIOBANK = 'MY_UKBIOBANK';
        DATASET_NAME = 'BATCH2_2021_DATA';
        GENO_ID = 'sel16875_rmrel';
end
RESULTS_ROOT = '../results/';
RESULTS_DIR = [RESULTS_ROOT, 'demoTorqueBrainAsymmetry/' DATASET_NAME '/'];
COMPUTED_ANGLES_PATH = [RESULTS_DIR, 'computedAngles.mat'];

THREADS = getenv('THREADS');
if(isempty(THREADS))
    THREADS = parcluster('local');
    THREADS = THREADS.NumWorkers;
else
    if ~isnumeric(THREADS)
        THREADS=str2double(THREADS);
    end
end
%% SETTING UP COMPUTATION POWER
try
    parpool(THREADS);
catch
end

COMPONENT=getenv('SYM_COMPONENT');
if isempty(COMPONENT)
    COMPONENT = DEFAULT_COMPONENT;
end

SELECTED_VARIANCE_THRESHOLD = 80;

disp(['Location of data: ', DATA_DIR]);
disp(['Using dataset:', num2str(DATASET_INDEX)])
disp(['Location of results: ', RESULTS_DIR]);

covGenoPath = [DATA_DIR, 'IMAGEN/BRAIN/' UKBIOBANK '/COVARIATES/COVDATAINLIERS.mat'];


if ~isfile(COMPUTED_ANGLES_PATH)

    in = load([DATA_DIR, 'HumanConnectomeProject/SubcorticalMask_HCP.mat']);

    phenopath = [DATA_DIR, 'IMAGEN/BRAIN/' UKBIOBANK '/PHENOTYPES/'];
    genopath = [DATA_DIR, 'IMAGEN/BRAIN/' UKBIOBANK '/GENOTYPES/'];

    Regions = {'LH' 'RH'};
    nR = length(Regions);

    DATA = cell(1,2);

    for r=1:nR
        regphenopath = [phenopath Regions{r} '/'];
        disp(['PROCESSING BRAIN REGION: ' Regions{r}]);
        DATA{r} = load([regphenopath DATASET_NAME '.mat']);
    end
    LH = DATA{1}.Region.AlignedShapes;
    RH = DATA{2}.Region.AlignedShapes;
    phenoIID = DATA{1}.Region.IID(1:size(LH,3));
    Region =  DATA{1}.Region;
    brainSurface = load([phenopath 'LH/RENDERMATERIAL.mat']);
    refTemplate = brainSurface.RefScan;

    % GPA
    [preprocTemplate, preprocLH, preprocRH, preprocPhenoIID, ~] = preprocessSymmetry(refTemplate, LH, RH, phenoIID, REDUCTION_RATE, 1, 3);
    %% Rerun once to allign to the template
    preprocRH(:,1,:) = -preprocRH(:,1,:);
    [preprocTemplate, preprocLH, preprocRH, preprocPhenoIID, preprocLandmarksIndices] = preprocessSymmetry(preprocTemplate, preprocLH, preprocRH, preprocPhenoIID, 1, 1, 1);

    %% Get torque
    close all
    angles = zeros(size(preprocLH,3),1);
    for ind=1:length(angles)
        indLH = preprocLH(:,:,ind);
        indRH = preprocRH(:,:,ind);
        indRH(:,1) = -indRH(:,1);
        [rhP, lhP, ang0] = getPoints(indRH, indLH);
        angles(ind) = ang0;
        if angles(ind) > 1.1 * pi/2 || angles(ind) < 0.9 * pi/2
            showTorque(indRH, indLH, preprocTemplate)
        end
    end
    avgLH = mean(preprocLH,3);
    avgRH = mean(preprocRH,3);
    avgRH(:,1) = -avgRH(:,1);
    save(COMPUTED_ANGLES_PATH, "angles", "preprocTemplate", 'avgRH', 'avgLH')
else
    load(COMPUTED_ANGLES_PATH)
end
%%
avgAngle = mean(angles);
fig = showTorque(avgRH, avgLH, preprocTemplate, avgAngle);

ang = 90 - rad2deg(avgAngle);
saveas(fig, sprintf('../results/torque_%.3f.svg',ang))

function fig =  showTorque(tRH, tLH, preprocTemplate, ang)
z_c = repmat(mean([tRH(:,3); tLH(:,3)]),size(tRH,1),1);

templateLH = clone(preprocTemplate);
templateRH = clone(preprocTemplate);
templateLH.Vertices = tLH;
templateRH.Vertices = tRH;

tRH = tRH(:,1:2);
tLH = tLH(:,1:2);
v1=templateLH.viewer;
templateRH.viewer(v1);
v1.BackgroundColor = 'white';
v1.SceneLightVisible = 1;
fig = v1.Figure;
if nargin < 4

    [rhP, lhP, ang0] = getPoints(tRH, tLH);
    ang = ang0;
    plot3(rhP(:, 1), rhP(:, 2), z_c(1:2), 'b', "LineWidth",4);
    plot3(lhP(:, 1), lhP(:, 2), z_c(1:2), 'b', "LineWidth",4);
end

ymin = min([tRH(:,2); tLH(:,2)]);
ymax = max([tRH(:,2); tLH(:,2)]);
ymid = (ymin + ymax)/2;
L = ymax - ymin;
x1 = 0 - (L/2 * cos(ang));
x2 = 0 + (L/2*cos(ang));
y1=ymid - (L/2*sin(ang));
y2=ymid+(L/2*sin(ang));



plot3(linspace(x1,x2, 1000),linspace(y1,y2, 1000),z_c(1:1000), 'r', "LineWidth",4);
plot3([0, 0], [y1, y2],z_c(1:2), 'k',"LineWidth",2)

end
%%
function [rhP, lhP, ang] = getPoints(RH, LH)
rLim = 1 / 3 * mean(RH(:,1));
lLim = 1 / 3 * mean(LH(:,1));
th = 0.01;
RH = RH(RH(:,1) < rLim, 1:2);
LH = LH(LH(:,1) > lLim, 1:2);
xyRH_hull = convhull(double(RH));
xyLH_hull = convhull(double(LH));

xyRH_hull_points = RH(xyRH_hull,1:2);
xyLH_hull_points = LH(xyLH_hull, 1:2);
rhPoints = xyRH_hull_points(:,1:2);
lhPoints = xyLH_hull_points(:,1:2);

rhDists = (rhPoints(2:end,1)-rhPoints(1:end-1,1)).^2 + (rhPoints(2:end,2)-rhPoints(1:end-1,2)).^2;
rhDists(rhPoints(2:end,1) > rLim - th) = 0;
lhDists = (lhPoints(2:end,1)-lhPoints(1:end-1,1)).^2 + (lhPoints(2:end,2)-lhPoints(1:end-1,2)).^2;
lhDists(lhPoints(2:end,1) < lLim + th) = 0;
[~, rhInd] = max(rhDists);
[~, lhInd] = max(lhDists);
rhP = rhPoints(rhInd:rhInd+1,:);
[~,a] = sort(rhP(:,2));
rhP = rhP(a,:);
lhP = lhPoints(lhInd:lhInd+1,:);
[~,b] = sort(lhP(:,2));
lhP = lhP(b,:);
rhAngleVec = [diff(rhP(:,1)), diff(rhP(:,2))];
lhAngleVec = [diff(lhP(:,1)), diff(lhP(:,2))];
avgVec = rhAngleVec/norm(rhAngleVec) + lhAngleVec/norm(lhAngleVec);
ang = angle(avgVec(1) + 1i * avgVec(2));
end