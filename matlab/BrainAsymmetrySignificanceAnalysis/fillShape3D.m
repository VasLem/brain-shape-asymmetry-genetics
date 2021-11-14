function reconstructed = fillShape3D(reconstructed, adjacency)
%fillShape3D in place fill Shape3D extracted Values, based on original
%adjacency matrix. Based on iterative smoothing (averaging diffusion)
% Adjacency is a matrix NxN (can be sparse)
% Reconstructed input has dimensionality (N,value_length) and
% downsampledIndices refer to the indices to which downsampledValues
% correspond to, will be inferred from the rest if not given
d2=size(reconstructed,2);
pdim = size(reconstructed,3);
for l = 1:50
    for c= 1:d2
    for k=1:pdim        
        if isnan(reconstructed(1, c,k))
           f = squeeze(reconstructed(:, c, adjacency(:,k) & ~ squeeze(isnan(reconstructed(1,c,:)))));
           if ~isempty(f)
               reconstructed(:, c, k) =  mean(f,2);
           end
        end
    end
    if sum(isnan(reconstructed))==0
        break
    end
    end
end
end