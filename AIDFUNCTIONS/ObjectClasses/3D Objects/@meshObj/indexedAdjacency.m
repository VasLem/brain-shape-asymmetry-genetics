function IA = indexedAdjacency(obj)
% generates Adjancency matrix with distances between connected points
         IA = obj.Adjacency;
         [i,j] = find(IA);
         ind = sub2ind(size(IA),i,j);
         IA(ind) = obj.IndexedColor;
         % make sure that all edges are symmetric
         %IA = max(IA,IA');
         IA = (IA+IA')/2;
end