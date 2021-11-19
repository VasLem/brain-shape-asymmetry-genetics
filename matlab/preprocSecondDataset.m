close all; clear
%% Phenotype
THREADS = 12;
try
    parpool('local',THREADS);
    %parpool('local',8);
catch
end
GPA_N = 3;
DATA_DIR = '../SAMPLE_DATA/';
W_OUT_PHENO_DIR = [DATA_DIR, 'IMAGEN/BRAIN/MY_UKBIOBANK/PHENOTYPES/'];
PHENO_DIR = '/usr/local/micapollo01/IMAGEN_DATA/UKbiobank/CIFTIFY/Batch2_2021/subjects/';
subDirs = dir(PHENO_DIR);
mask = load('/usr/local/micapollo01/MIC/DATA/STUDENTS/vlemon0/code/SAMPLE_DATA/IMAGEN/BRAIN/HumanConnectomeProject/SubcorticalMask_HCP.mat').index;
sides = ["L", "R"];
sidesRet = zeros(length(mask), 3, length(subDirs));
Regions = {'LH' 'RH'};
for sideInd = 1:2
    disp(['Processing region ' Regions{sideInd} '..'])
    disp("Reading using gifti..")
    ppb = ParforProgressbar(length(subDirs));
    parfor subDirInd=1:length(subDirs)
        if strcmp(subDirs(subDirInd).name(1), '.'), continue; end
        subject = subDirs(subDirInd).name;
        iid = split(subject,'_');
        iids{subDirInd} = iid(1);
        path = sprintf('%s/%s/MNINonLinear/fsaverage_LR32k/%s.%s.midthickness.32k_fs_LR.surf.gii',PHENO_DIR, subject, subject, sides(sideInd));
        ret = readSurface(path, mask, 0);
        sidesRet(:, :, subDirInd) = ret;
        ppb.increment();
    end
    delete(ppb);
    disp("Applying GPA..")
    subject = subDirs(10).name;
    path = sprintf('%s/%s/MNINonLinear/fsaverage_LR32k/%s.%s.midthickness.32k_fs_LR.surf.gii',PHENO_DIR, subject, subject, sides(sideInd));
    template = readSurface(path, mask, 1);
    
    [Region.AlignedShapes,~,~] = GeneralizedProcrustesAnalysis(sidesRet, template, GPA_N,true,false,true,false);
    template.Vertices = mean(Region.AlignedShapes, 3);
    Region.AvgShape = template;
    Region.IID = iids';
    disp("Saving..")
    outdir = [W_OUT_PHENO_DIR Regions{sideInd} '/'];
    if ~isfolder(outdir), mkdir(outdir); end
    outpath = [outdir 'BATCH2_2021_DATA'];
    save(outpath,'Region','-v7.3');
end

%%Genotype




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
