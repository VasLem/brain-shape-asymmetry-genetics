function [in,prob] = inClusters(clusters,uv,th)
% returns the index of the cluster a point belongs to 
% and 0 if the probabilty is to low according to th 
% to belong to any of the clusters
     if size(uv,1) == 2;
        uv = uv';
     end
     nrC = length(clusters);
     y = zeros(nrC,size(uv,1));
     for k=1:1:nrC
         y(k,:) = mvnpdf(uv,clusters(k).mu',clusters(k).cov)./clusters(k).max_prob;         
     end
     if nrC>1
        [maxval,in] = max(y);
        in(find(maxval<th)) = 0; %#ok<FNDSB>
        prob = maxval;
     else
         in = ones(size(y));
         in(find(y<th)) = 0; %#ok<FNDSB>
         prob = y;
     end
        
end
