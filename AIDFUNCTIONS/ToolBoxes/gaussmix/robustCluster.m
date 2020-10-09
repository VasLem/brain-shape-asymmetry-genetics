function [cluster,uv,index] = robustCluster(uv,index)
         nrE = length(index);
         subindex = [];
         counter = 0;
         while ~(nrE==length(subindex))&&counter<20
             counter = counter +1;
             nrE = length(index);
             %disp(num2str(counter));
             mtrs = GaussianMixture(uv(:,index)',1,1,false);
             cluster = getClusterInfo(mtrs,uv(:,index));
             y = mvnpdf(uv(:,index)',[0.5 0.5],[0.5 0;0 0.5]);
             subindex = find(cluster.y*cluster.max_prob>y);
             index = index(subindex);
         end
         mtrs = GaussianMixture(uv(:,index)',1,1,false);
         cluster = getClusterInfo(mtrs,uv(:,index));
         uv = uv(:,index);         
end
