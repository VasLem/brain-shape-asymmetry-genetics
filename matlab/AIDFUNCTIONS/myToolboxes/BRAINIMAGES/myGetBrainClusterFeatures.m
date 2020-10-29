function [Size,Shape] = myGetBrainClusterFeatures(data)

    [nObs,nVoxels] = size(data);
    Size = sum(data,2)./nVoxels;





end