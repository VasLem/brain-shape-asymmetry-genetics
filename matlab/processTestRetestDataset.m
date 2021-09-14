addpath /home/vlemon0/home/gifti
addpath('AIDFUNCTIONS');
close all

in = load('/usr/local/micapollo01/MIC/DATA/STUDENTS/vlemon0/code/SAMPLE_DATA/IMAGEN/BRAIN/HumanConnectomeProject/SubcorticalMask_HCP.mat');
for subNum=1:20
    for session=1:2
        dataset.LH(:, :, subNum,session) = readSurface(num2str(subject, '%02.f'), 'L',  num2str(session), in.index).Vertices;
        dataset.RH(:, :, subNum,session) = readSurface(num2str(subject, '%02.f'), 'R', num2str(session), in.index).Vertices;
    end
end
%%

shape = size(dataset.LH);
%%
brainSurface = load('../SAMPLE_DATA/IMAGEN/BRAIN/UKBIOBANK/PHENOTYPES/LH/RENDERMATERIAL.mat');
pLH = permute(reshape(permute(dataset.LH, [3,4, 1,2]), [shape(3) * shape(4), shape(1:2)]), [2, 3, 1]);
pRH = permute(reshape(permute(dataset.RH, [3,4, 1,2]), [shape(3) * shape(4), shape(1:2)]), [2, 3, 1]);
[~, alignedPLH, alignedPRH, ~, ~] = preprocessSymmetry(brainSurface.RefScan, pLH, pRH, (1:size(pLH,3)));
alignedLH = permute(reshape(permute(alignedPLH, [3, 1, 2]), [shape(3), shape(4), shape(1), shape(2)]), [3,4,1,2]);
alignedRH = permute(reshape(permute(alignedPRH, [3, 1, 2]), [shape(3), shape(4), shape(1), shape(2)]), [3,4,1,2]);
%%
outLH = mean(std(alignedLH,0,4),3);
outRH = mean(std(alignedRH,0,4),3);
%The following is the reason why we decided to keep both variances and not
%compute the mean out of them. For the test retest dataset, technical
%measurement errors do not seem to be characterized by symmetry on the
%midsaggital plane
fig = figure;
hist(outLH - outRH);
title('Difference between measurement variances of contralateral and symmetrical, as of the midsaggital plane, landmarks')
variance.LH = outLH;
variance.RH = outRH;
%%
save("../SAMPLE_DATA/test_retest_dataset.mat", "dataset", "-v7");
save("../results/test_retest_information.mat", "variance","-v7");
%%

function ret = readSurface(subject, side, session, mask)
    direc = '/usr/local/micapollo01/IMAGEN_DATA/OTHER/IMAGEDATA/BRAINDATA/ReplicationStudy_Radwan/Ciftify/';
    subSesID = ['sub-S0' subject '_ses-0' session '_brain'];
    subSesDir = [direc subSesID '/MNINonLinear/fsaverage_LR32k/'];
    path = [subSesDir  subSesID '.' side '.midthickness.32k_fs_LR.surf.gii'];
    g = gifti(path);
    ret = shape3D;
    ret.Vertices = g.vertices;
    ret.Faces = g.faces;
    ret = crop(ret, 'VertexIndex', mask);
end