
function dataset = loadTestRetestDataset()
try
    load("../SAMPLE_DATA/test_retest_dataset.mat", "dataset")
catch
    addpath /home/vlemon0/home/gifti
    in = load('/usr/local/micapollo01/MIC/DATA/STUDENTS/vlemon0/code/SAMPLE_DATA/IMAGEN/BRAIN/HumanConnectomeProject/SubcorticalMask_HCP.mat');
    for subNum=1:20
        for session=1:2
            dataset.LH(:, :, subNum,session) = readSurface(num2str(subNum, '%02.f'), 'L',  num2str(session), in.index).Vertices;
            dataset.RH(:, :, subNum,session) = readSurface(num2str(subNum, '%02.f'), 'R', num2str(session), in.index).Vertices;
        end
        
    end
    
    save("../SAMPLE_DATA/test_retest_dataset.mat", "dataset", "-v7");
end
end


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