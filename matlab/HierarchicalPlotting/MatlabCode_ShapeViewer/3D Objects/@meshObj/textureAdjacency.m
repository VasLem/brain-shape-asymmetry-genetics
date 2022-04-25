function TA = textureAdjacency(obj)
% generates Adjancency matrix with distances between connected points
         TA = obj.Adjacency;
         [i,j] = find(TA);
         ind = sub2ind(size(TA),i,j);
         distances = sqrt(sum((obj.TextureColor(:,i)-obj.TextureColor(:,j)).^2));
         TA(ind) = distances;
         % make sure that all edges are symmetric
         TA = (TA+TA')/2;
end