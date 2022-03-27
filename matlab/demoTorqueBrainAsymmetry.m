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
        xyRH_hull = convhull(double(indRH(:,1:2)));
        xyLH_hull = convhull(double(indLH(:,1:2)));
        xyRH_hull_points = indRH(xyRH_hull,1:2);
        xyLH_hull_points = indLH(xyLH_hull, 1:2);
        rhPoints = xyRH_hull_points(xyRH_hull_points(:,1)<0.1,:);
        lhPoints = xyLH_hull_points(xyLH_hull_points(:,1)>-0.1,:);

        rhDists = (rhPoints(2:end,1)-rhPoints(1:end-1,1)).^2 + (rhPoints(2:end,2)-rhPoints(1:end-1,2)).^2;
        lhDists = (lhPoints(2:end,1)-lhPoints(1:end-1,1)).^2 + (lhPoints(2:end,2)-lhPoints(1:end-1,2)).^2;
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
        p = [0; 0; 1];
        %     angdiff = vecangle360(rhAngleVec,lhAngleVec,p);
        %     avgAngle = angle(rhAngleVec(1) + 1i * rhAngleVec(2))+ angdiff / 2;
        %     angles(ind) = avgAngle;
        avgVec = rhAngleVec/norm(rhAngleVec) + lhAngleVec/norm(lhAngleVec);
        angles(ind) = angle(avgVec(1) + 1i * avgVec(2));
    end
    avgLH = mean(preprocLH,3);
    avgRH = mean(preprocRH,3);
    avgRH(:,1) = -avgRH(:,1);
    save(COMPUTED_ANGLES_PATH, "angles", "preprocTemplate", 'avgRH', 'avgLH')
else
    load(COMPUTED_ANGLES_PATH)
end



templateLH = clone(preprocTemplate);
templateRH = clone(preprocTemplate);

templateLH.Vertices = avgLH;
templateRH.Vertices = avgRH;
v1=templateLH.viewer;
templateRH.viewer(v1);
v1.BackgroundColor = 'white';
v1.SceneLightVisible = 1;

xyRH_hull = convhull(double(avgRH(:,1:2)));
xyLH_hull = convhull(double(avgLH(:,1:2)));
%%


% plot3(avgRH(xyRH_hull,1), avgRH(xyRH_hull,2),z_c(xyRH_hull),'r')
% plot3(avgLH(xyLH_hull,1), avgLH(xyLH_hull,2),z_c(xyLH_hull),'r')
tRH = avgRH;
tLH = avgLH;
% tRH = indRH;
% tLH = indLH;
xyRH_hull = convhull(double(tRH(:,1:2)));
xyLH_hull = convhull(double(tLH(:,1:2)));
xyRH_hull_points = tRH(xyRH_hull,:);
xyLH_hull_points = tLH(xyLH_hull, :);
rhPoints = xyRH_hull_points(xyRH_hull_points(:,1)<0.1,:);
lhPoints = xyLH_hull_points(xyLH_hull_points(:,1)>-0.1,:);
rhDists = (rhPoints(2:end,1)-rhPoints(1:end-1,1)).^2 + (rhPoints(2:end,2)-rhPoints(1:end-1,2)).^2;
lhDists = (lhPoints(2:end,1)-lhPoints(1:end-1,1)).^2 + (lhPoints(2:end,2)-lhPoints(1:end-1,2)).^2;
[~, rhInd] = max(rhDists);
[~, lhInd] = max(lhDists);
z_c = repmat(mean([tRH(:,3); tLH(:,3)]),size(tRH,1),1);
scatter3(rhPoints(:,1), rhPoints(:,2), z_c(1:length(rhPoints)), 'b', '.');
scatter3(lhPoints(:,1), lhPoints(:,2), z_c(1:length(lhPoints)), 'b', '.');

plot3(rhPoints(rhInd:rhInd+1,1), rhPoints(rhInd:rhInd+1,2), z_c(1:2), 'b');
plot3(lhPoints(lhInd:lhInd+1,1), lhPoints(lhInd:lhInd+1,2), z_c(1:2), 'b');
avgAngle = mean(angles);

ymin = min([xyRH_hull_points(:,2); xyLH_hull_points(:,2)]);
ymax = max([xyRH_hull_points(:,2); xyLH_hull_points(:,2)]);
ymid = (ymin + ymax)/2;
L = ymax - ymin;
x1 = 0 - (L/2 * cos(avgAngle));
x2 = 0 + (L/2*cos(avgAngle));
y1=ymid - (L/2*sin(avgAngle));
y2=ymid+(L/2*sin(avgAngle));
plot3(linspace(x1,x2, 1000),linspace(y1,y2, 1000),z_c(1:1000), 'r', "LineWidth",4);
plot3([0, 0], [y1, y2],z_c(1:2), 'k',"LineWidth",2)
%%
ang = 90 + rad2deg(avgAngle);
saveas(v1.Figure, sprintf('../results/torque_%.3f.svg',ang))