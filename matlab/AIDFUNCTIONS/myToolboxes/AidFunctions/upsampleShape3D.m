
function reconstructed = upsampleShape3D(downsampledValues, adjacency, downsampledIndices)
%UPSAMPLESHAPE3D Upsample Shape3D extracted Values, based on original
%adjacency matrix. Based on iterative smoothing (averaging diffusion)
% Adjacency is a matrix NxN (can be sparse)
% Expected values have dimensionality (reducedN,value_length) and
% downsampledIndices refer to the indices to which downsampledValues
% correspond to, will be inferred from the rest if not given

pdim = size(adjacency,1);
if length(size(downsampledValues)) == 2
    downsampledValues = reshape(downsampledValues,[size(downsampledValues,1),1,size(downsampledValues,2)]);
end
odim = size(downsampledValues,3);
if nargin < 3, downsampledIndices = 1: ceil(pdim/odim) : pdim ; end
d1=size(downsampledValues, 1);
d2=size(downsampledValues,2);

reconstructed = zeros(d1,d2, pdim);
reconstructed(:) = NaN;
reconstructed(:, :, downsampledIndices) = downsampledValues;
reconstructed = fillShape3D(reconstructed, adjacency);
end


