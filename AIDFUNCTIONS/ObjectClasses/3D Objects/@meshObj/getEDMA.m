function M = getEDMA(obj)
         M = squareform(pdist(obj.Vertices','euclidean'));
end