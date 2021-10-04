function FA = functionAdjacency(obj,f)
% generates Adjancency matrix with distances between connected points
         FA = obj.Adjacency;
         [i,j] = find(FA);
         ind = sub2ind(size(FA),i,j);
         distances = sqrt(sum((f(i)-f(j)).^2));
         FA(ind) = distances;
         % make sure that all edges are symmetric
         %FA = max(FA,FA');
         FA = (FA+FA')/2;
end