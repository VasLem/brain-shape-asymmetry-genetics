function DA = vertexAdjacency(obj)
% generates Adjancency matrix with distances between connected points
         DA = obj.Adjacency;
         [i,j] = find(DA);
         ind = sub2ind(size(DA),i,j);
         distances = sqrt(sum((obj.Vertices(:,i)-obj.Vertices(:,j)).^2));
         DA(ind) = distances;
         % make sure that all edges are symmetric
         DA = (DA+DA')/2;
end