clc;clear all;
addpath(genpath('/usr/local/micapollo01/MIC/DATA/STAFF/myuan0/tmp/BrainG2L/SHARED/AIDFUNCTIONS/'));
datapath = '/IMAGEN/BRAIN/UKBIOBANK/CIFTIFY/';
outpath = '/usr/local/micapollo01/MIC/DATA/STAFF/myuan0/tmp/BrainG2L/01_Pheno/Data/Stage00Data/';
 
load('/IMAGEN/BRAIN/UKBIOBANK/COVARIATES/COVDATAINLIERS.mat');
% histcounts(COV.DATA(:,21))
in = load('/IMAGEN/BRAIN/UKBIOBANK/CIFTIFY_SAMPLE/SubcorticalMask_HCP.mat');
MASK = in.index;
nVertices = length(MASK);
%%
IID = COV.IID;% These IDs are selected, starting from images (20407, no outliers and artifacts) to genetics (no outliers and relatives)
nI = length(IID);
Regions = {'R' 'L'};
nR = length(Regions);
for r=1:1:nR
    Region.Name = Regions{r};
    disp(['PROCESSING BRAIN REGION: ' Region.Name]);
    present = zeros(1,nI);
    THICKNESS = nan*zeros(nVertices,1,nI,'single');  
    [path,ID] = setupParForProgress(nI);
    parfor i=1:nI
        imgpath = [datapath IID{i} '_brain/MNINonLinear/fsaverage_LR32k/'];
        
        imgpath_piral = [imgpath IID{i} '_brain.' Region.Name '.pial.32k_fs_LR.surf.gii'];
        g_piral = gifti(imgpath_piral);
        obj2 = shape3D;
        obj2.Vertices = double(g_piral.vertices);
        obj2.Faces = double(g_piral.faces);
        crop(obj2,'VertexIndex',MASK);
        
        imgpath_white = [imgpath IID{i} '_brain.' Region.Name '.white.32k_fs_LR.surf.gii'];
        g_white = gifti(imgpath_white);
        obj3 = shape3D;
        obj3.Vertices = double(g_white.vertices);
        obj3.Faces = double(g_white.faces);
        crop(obj3,'VertexIndex',MASK);
        
        present(i) = 1;
        thickness = sqrt((obj2.Vertices(:,1)-obj3.Vertices(:,1)).^2+...
            (obj2.Vertices(:,2)-obj3.Vertices(:,2)).^2+...
            (obj2.Vertices(:,3)-obj3.Vertices(:,3)).^2);       
        THICKNESS(:,:,i) = single(thickness);
        parfor_progress;
    end
    closeParForProgress(path,ID);
   
    index = find(present);
    Region.IID = IID(index); 
    THICKNESS = THICKNESS(:,:,index);
    save([outpath 'Stage00Thickness_' Regions{r} 'H'],'THICKNESS','-v7.3');
end

STAGE00THICKNESS_L = load([outpath '/Stage00Thickness_LH.mat']);
Thickness_L = STAGE00THICKNESS_L.THICKNESS;
STAGE00THICKNESS_R = load([outpath '/Stage00Thickness_RH.mat']);
Thickness_R = STAGE00THICKNESS_R.THICKNESS;
% average of left/right brain
Stage00ThicknessSym = (Thickness_L + Thickness_R)/2;
save([outpath 'Stage00ThicknessSym'],'Stage00ThicknessSym','-v7.3');


