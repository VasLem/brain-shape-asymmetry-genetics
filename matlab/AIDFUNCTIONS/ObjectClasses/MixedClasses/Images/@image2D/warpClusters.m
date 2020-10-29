function warpedclusters = warpClusters(obj,clusters,startUV,newUV)
         warning off; %#ok<WNOFF>
         %newUV = obj2.UV; obj = obj1;
         nrC = length(clusters);
         in = inClusters(clusters,startUV,0.01);        
         for i=1:1:nrC
             Cindex = find(in==i);
             [cluster,newCUV,Cindex] = robustCluster(newUV,Cindex);
             ProbRBF = getProbRBF(newCUV,cluster.y);
             warped = warpClusterImage(obj,startUV(:,Cindex),newCUV);             
             eval = fastrbf_pointeval(ProbRBF,warped.UV,'messages',0);
             warpedclusters(i).Warped = warped;clear warped; %#ok<AGROW>
             warpedclusters(i).ProbRBF = ProbRBF;clear ProbRBF prob; %#ok<AGROW>
             warpedclusters(i).Prob = eval.Value;clear eval; %#ok<AGROW>
         end
end

% function ProbRBF = getProbRBF(UV, prob)
%          % UV = newCUV;
%          pl.Location = UV;
%          pl.Value = prob';
%          %pl = fastrbf_unique(pl,'messages',0);
%          pl = uniquePL(pl);
%          ProbRBF = fastrbf_fit(pl,0.01,'reduce','messages',0);
% end

% function [in,prob] = inClusters(clusters,uv,th)
% % returns the index of the cluster a point belongs to 
% % and 0 if the probabilty is to low according to th 
% % to belong to any of the clusters
%      if size(uv,1) == 2;
%         uv = uv';
%      end
%      nrC = length(clusters);
%      y = zeros(nrC,size(uv,1));
%      for k=1:1:nrC
%          y(k,:) = mvnpdf(uv,clusters(k).mu',clusters(k).cov)./clusters(k).max_prob;         
%      end
%      if nrC>1
%         [maxval,in] = max(y);
%         in(find(maxval<th)) = 0; %#ok<FNDSB>
%         prob = maxval;
%      else
%          in = ones(size(y));
%          in(find(y<th)) = 0; %#ok<FNDSB>
%          prob = y;
%      end
%         
% end


% function [cluster,uv,index] = robustCluster(uv,index)
%          nrE = length(index);
%          subindex = 1;
%          counter = 0;
%          while ~(nrE==length(subindex))&&counter<10
%              counter = counter +1;
%              nrE = length(index);
%              %disp(num2str(counter));
%              mtrs = GaussianMixture(uv(:,index)',1,1,false);
%              cluster = getClusterInfo(mtrs,uv(:,index));
%              y = mvnpdf(uv(:,index)',[0.5 0.5],[0.5 0;0 0.5]);
%              subindex = find(cluster.y*cluster.max_prob>y);
%              index = index(subindex);
%          end
%          uv = uv(:,index);         
% end

