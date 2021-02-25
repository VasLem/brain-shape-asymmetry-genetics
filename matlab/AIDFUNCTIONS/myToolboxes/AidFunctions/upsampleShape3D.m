
function reconstructed = upsampleShape3D(downsampledValues, adjacency, downsampledIndices)
%UPSAMPLESHAPE3D Upsample Shape3D extracted Values, based on original
%adjacency matrix. Based on iterative smoothing (averaging diffusion)
% Adjacency is a matrix NxN (can be sparse)
% Expected values have dimensionality (reducedN,value_length) and
% downsampledIndices refer to the indices to which downsampledValues
% correspond to, will be inferred from the rest if not given

if nargin < 3, downsampledIndices = 1: ceil(size(adjacency,1) / size(downsampledValues,1)) :size(adjacency,1) ; end
pdim = size(adjacency,1);
reconstructed = zeros(size(adjacency,1),size(downsampledValues,2));
reconstructed(:) = NaN;
reconstructed(downsampledIndices, :) =downsampledValues;
for l = 1:50
    for k=1:pdim        
        if isnan(reconstructed(k,1))
           f = reconstructed(adjacency(:,k) & ~ isnan(reconstructed(:,1)),:);
           if ~isempty(f)
               reconstructed(k,:) =  mean(f,1);
           end
        end
    end
    if sum(isnan(reconstructed))==0
        break
    end
end