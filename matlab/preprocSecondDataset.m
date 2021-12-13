close all;clear;
restoredefaultpath;
addpath /home/vlemon0/home/gifti
addpath(genpath('AIDFUNCTIONS'));
DATA_DIR = '../SAMPLE_DATA/';
THREADS =8;
try
    parpool('local',THREADS);
    %parpool('local',8);
catch
end
addpath(genpath('/opt/SNPLIB/'));
rmpath('/opt/SNPLIB/mexfiles/');% to remove the functions that are causing matlab to crash
addpath(genpath('SNPLIB-master/mexfiles/'))% where I stored the re-mexed files
%% Phenotype
GPA_N = 3;
W_OUT_PHENO_DIR = [DATA_DIR, 'IMAGEN/BRAIN/MY_UKBIOBANK/PHENOTYPES/'];
PHENO_DIR = '/usr/local/micapollo01/IMAGEN_DATA/UKbiobank/CIFTIFY/Batch2_2021/subjects/';
subDirs = dir(PHENO_DIR);
mask = load('/usr/local/micapollo01/MIC/DATA/STUDENTS/vlemon0/code/SAMPLE_DATA/IMAGEN/BRAIN/HumanConnectomeProject/SubcorticalMask_HCP.mat').index;
sides = ["L", "R"];
Regions = {'LH' 'RH'};
for sideInd = 1:2
    disp(['Processing region ' Regions{sideInd} '..'])
    disp("Reading using gifti..")
    ppb = ParforProgressbar(length(subDirs));
    sidesRet = zeros(length(mask), 3, length(subDirs));
    centroidSizes = zeros(length(subDirs), 1);
    iids = cell(length(subDirs), 1);
    parfor subDirInd=1:length(subDirs)
        if strcmp(subDirs(subDirInd).name(1), '.'), continue; end
        subject = subDirs(subDirInd).name;
        iid = split(subject,'_');
        iids{subDirInd} = iid(1);
        path = sprintf('%s/%s/MNINonLinear/fsaverage_LR32k/%s.%s.midthickness.32k_fs_LR.surf.gii',PHENO_DIR, subject, subject, sides(sideInd));
        ret = readSurface(path, mask, 0);
        sidesRet(:, :, subDirInd) = ret;
        centroidSizes(subDirInd) = computeCentroidSize(ret);
        ppb.increment();
    end
    check = ~cellfun(@isempty, iids);
    iids = iids(check);
    iids = cellfun(@char, iids, 'UniformOutput', false);

    sidesRet = sidesRet(:, :, check);
    centroidSizes = centroidSizes(check);
    delete(ppb);
    disp("Applying GPA..")
    subject = subDirs(10).name;
    path = sprintf('%s/%s/MNINonLinear/fsaverage_LR32k/%s.%s.midthickness.32k_fs_LR.surf.gii',PHENO_DIR, subject, subject, sides(sideInd));
    template = readSurface(path, mask, 1);
    
    [Region.AlignedShapes,~,~] = GeneralizedProcrustesAnalysis(sidesRet, template, GPA_N,true,false,true,false);
    template.Vertices = mean(Region.AlignedShapes, 3);
    Region.AvgShape = template;
    Region.Name =  Regions{sideInd};
    Region.CentroidSizes = centroidSizes;
    if sideInd == 1
        o = iids;
    else
        assert(all(str2double(o) == str2double(iids)));
    end
    Region.IID = iids;
    disp("Saving..")
    outdir = [W_OUT_PHENO_DIR Regions{sideInd} '/'];
    if ~isfolder(outdir), mkdir(outdir); end
    outpath = [outdir 'BATCH2_2021_DATA'];
    save(outpath,'Region','-v7.3');
    save([W_OUT_PHENO_DIR, 'BATCH2_2021_DATA_IID'], 'iids', '-v7.3');
end

%% Genotype
disp("Moving on to Genotype...")
W_OUT_GENO_DIR = [DATA_DIR, '../SAMPLE_DATA/IMAGEN/BRAIN/MY_UKBIOBANK/GENOTYPES/PLINK/'];
GENO_DIR = '/usr/local/micapollo01/IMAGEN_DATA/SHARED/sgoova5/UKB_batch2_genotypes/3_RMREL/';
GENO_ID = 'ukb_img_maf0.01_geno0.5_hwe1e-6_sel16875_rmrel';
GENO_FILE_SUFFIX = '_ALLchr';
PLINK_PATH = '../bash/genomics/plink';

bfile = sprintf('%s%s%s',GENO_DIR, GENO_ID, GENO_FILE_SUFFIX );
if ~isfolder(W_OUT_GENO_DIR), mkdir(W_OUT_GENO_DIR); end
disp("Making PCA components file..")
cmd = sprintf('%s --noweb --bfile %s  --memory 32000 --threads 16 --pca --out %s%s%s',PLINK_PATH,bfile , W_OUT_GENO_DIR, GENO_ID,'_pca');
system(cmd);

% geno(geno==-1) = 255;
% geno = uint8(geno);
%%
AA= pca(geno, 'Centered', true, 'NumComponents', 20);
FID =samples.IID;
IID = samples.IID;
save('../SAMPLE_DATA/IMAGEN/BRAIN/MY_UKBIOBANK/GENOTYPES/UKB_EUR_16875_PCs_2.mat', 'AA', 'FID', 'IID', '-v7.3');

%%
disp("Splitting chromosomes..")
parfor chr=1:22
    cmd = sprintf('%s --noweb --bfile %s --chr %s --make-bed --out %s%s_chr%s',PLINK_PATH,bfile , num2str(chr), W_OUT_GENO_DIR, GENO_ID,  num2str(chr));
    system(cmd);
end



function ret = readSurface(path, mask, retShape)
g = gifti(path);
if ~retShape
ret  = g.vertices(mask, :);
else
ret = shape3D;
ret.Vertices = g.vertices;
ret.Faces = g.faces;
ret = crop(ret, 'VertexIndex', mask);
end
end

function out = computeCentroidSize(verts)
    nVertices = size(verts,1);
    centroid =  mean(verts,1); 
    Differences = repmat(centroid,nVertices,1)-nVertices;
    Distances = sqrt(sum(Differences.^2,2));
    out = sqrt(sum(Distances.^2)/nVertices);
end