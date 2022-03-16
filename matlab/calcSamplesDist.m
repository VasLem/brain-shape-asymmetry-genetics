
close all;clear;
if ~isdeployed
    restoredefaultpath;
    addpath(genpath('AIDFUNCTIONS'));
    addpath(genpath('FUNCTIONS'));
    addpath(genpath('.'));
end
NUM_LEVELS = 4;

DEFAULT_DATASET_INDEX = 1;
REDUCTION_RATE = 0.1;
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
THREADS = parcluster('local');
THREADS = THREADS.NumWorkers;

%% SETTING UP COMPUTATION POWER
try
    parpool(THREADS);
catch
end


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
brainSurface = load([regphenopath 'RENDERMATERIAL.mat']);
refTemplate = brainSurface.RefScan;

% GPA
[preprocTemplate, preprocLH, preprocRH, preprocPhenoIID, preprocLandmarksIndices] = preprocessSymmetry(refTemplate, LH, RH, phenoIID, REDUCTION_RATE, 1, 3);
%%
distmats = cell(2,1);
for c=1:2
    if c == 1
        mat = preprocLH;
        disp('Processing left');
    else
        mat = preprocRH;
        disp('Processing right');
    end
    ret = cell(size(mat,3),1);
    p =  ParforProgressbar(size(mat,3));
    parfor i=1:size(mat,3)
        diff = mat(:,:,i+1:end) - mat(:, :, i);
        diff = mean(sqrt(diff(:,1,:).^2 + diff(:,2,:).^2 + diff(:,3,:).^2));
        r = zeros(1, size(mat,3));
        r(i+1:end) = diff;
        ret{i} = r;
        p.increment();
    end
    delete(p); 
    distmats{c} = cell2mat(ret);
    distmats{c} = distmats{c} + distmats{c}';
end

