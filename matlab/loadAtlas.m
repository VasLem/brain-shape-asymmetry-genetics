function atlas = loadAtlas(atlasName, side)
    if nargin < 2, side = 'L'; end
    load('../SAMPLE_DATA/IMAGEN/BRAIN/HumanConnectomeProject/SubcorticalMask_HCP.mat', 'index');
    mask = index;
    atlas = load(['../SAMPLE_DATA/atlasses/', side , atlasName, '_Average.mat'], atlasName);
    atlas = atlas.(atlasName);
    atlas.index = atlas.index(mask);
end