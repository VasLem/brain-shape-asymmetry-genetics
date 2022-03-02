clear
close all
DATA_DIR = getenv('DATA_ROOT');
if(isempty(DATA_DIR))
    DATA_DIR = '../SAMPLE_DATA/';
end
DATASET_INDEX = getenv('DATASET_INDEX');
if (isempty(DATASET_INDEX))
    DATASET_INDEX = 1;
else
    disp(DATASET_INDEX)
    if ~isnumeric(DATASET_INDEX)
        DATASET_INDEX = str2double(DATASET_INDEX);
    end
end

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
in = load([DATA_DIR, '/HumanConnectomeProject/SubcorticalMask_HCP.mat']);

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
[preprocTemplate, preprocLH, preprocRH, preprocPhenoIID, preprocLandmarksIndices] = preprocessSymmetry(refTemplate, LH, RH, phenoIID, 1, 1, 3);

%% Get torque
avgLH = mean(LH,3);
avgRH = mean(RH,3);
%%
templateLH = clone(preprocTemplate);
templateRH = clone(preprocTemplate);
%%
templateLH.Vertices = avgLH;
templateRH.Vertices = avgRH;
%%
v1=templateLH.viewer;
templateRH.viewer(v1);
v1.BackgroundColor = 'white';
v1.SceneLightVisible = 1;
xyRH_hull = convhull(double(avgRH(:,1:2)));
xyLH_hull = convhull(double(avgLH(:,1:2)));
z_c = repmat(mean([avgRH(:,3); avgLH(:,3)]),size(avgRH,1),1);

% plot3(avgRH(xyRH_hull,1), avgRH(xyRH_hull,2),z_c(xyRH_hull),'r')
% plot3(avgLH(xyLH_hull,1), avgLH(xyLH_hull,2),z_c(xyLH_hull),'r')

xyRH_hull_points = avgRH(xyRH_hull,:);
xyLH_hull_points = avgLH(xyLH_hull, :);
rhPoints = xyRH_hull_points(xyRH_hull_points(:,1)<0.08,:);
lhPoints = xyLH_hull_points(xyLH_hull_points(:,1)>-0.08,:);
scatter3(rhPoints(:,1), rhPoints(:,2), z_c(1:length(rhPoints)), 'b', '.');
scatter3(lhPoints(:,1), lhPoints(:,2), z_c(1:length(lhPoints)), 'b', '.');

rhDists = (rhPoints(2:end,1)-rhPoints(1:end-1,1)).^2 + (rhPoints(2:end,2)-rhPoints(1:end-1,2)).^2;
lhDists = (lhPoints(2:end,1)-lhPoints(1:end-1,1)).^2 + (lhPoints(2:end,2)-lhPoints(1:end-1,2)).^2;

[~, rhInd] = max(rhDists);
[~, lhInd] = max(lhDists);

plot3(rhPoints(rhInd:rhInd+1,1), rhPoints(rhInd:rhInd+1,2), z_c(1:2), 'b');
plot3(lhPoints(lhInd:lhInd+1,1), lhPoints(lhInd:lhInd+1,2), z_c(1:2), 'b');

rhAngle = atan(diff(rhPoints(rhInd:rhInd+1,2)) / diff(rhPoints(rhInd:rhInd+1,1)));
lhAngle = atan(diff(lhPoints(lhInd:lhInd+1,2)) / diff(lhPoints(lhInd:lhInd+1,1)));
avgAngle = mean([rhAngle, lhAngle]);

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