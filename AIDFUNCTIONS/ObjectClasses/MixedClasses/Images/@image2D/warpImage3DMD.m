function out = warpImage3DMD(obj,startUV,newUV,res)
         if nargin < 4
            res = [0.002 0.002];
         end
         if nargout == 1
            obj = clone(obj);
            out = obj;
         end
%          backupobj = obj;
%          obj = Target.TextureMap;
%          startUV = Target.UV;
%          newUV = UV;
%          res = [0.002 0.002];
         % Getting Clusters
         mtrs = GaussianMixture(startUV', 4, 2, false);
         clusters = getClusterInfo(mtrs(1),startUV);             
%          in = inClusters(clusters,startUV,0.01);
%          figure;fastrbf_view(startUV,in);
%          figure;fastrbf_view(newUV,in);
%          i=1;
%          Cindex = find(in==i);
%          [cluster,newCUV,Cindex] = robustCluster(newUV,Cindex);
%          figure;fastrbf_view(newCUV,cluster.y');
%          
%          uv = newUV;
%          index = Cindex;
%          nrE = length(index);
%          subindex = [];
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
                
         % warping clusters
         clusters = warpClusters(obj,clusters,startUV,newUV);
%          for i=1:1:2
%              figure;fastrbf_view(clusters(i).Warped.UV,clusters(i).Warped.Image);
%              figure;fastrbf_view(clusters(i).Warped.UV,clusters(i).Prob);
%          end
         clusters = clusterInfluence(clusters);
%          for i=1:1:2
%              figure;fastrbf_view(clusters(i).Warped.UV,clusters(i).Warped.Image);
%              figure;fastrbf_view(clusters(i).Warped.UV,clusters(i).Prob);
%          end
         Location = [];
         Value = [];
         for c=1:1:length(clusters)
             Location = [Location clusters(c).Warped.UV]; %#ok<AGROW>
             Value = [Value clusters(c).Warped.Image]; %#ok<AGROW>
         end        
         [X Y] = meshgrid(0:res(1):1,0:res(2):1);
         im = zeros(size(X,1),size(X,2),obj.Dim);
         warning('off','All');
         for i=1:1:obj.Dim
             im(:,:,i) = griddata(Location(1,:),Location(2,:),Value(i,:),X,Y,'cubic');
         end
%          figure;fastrbf_view(Location,Value);
         clear warped;
         warning('on','All');
         obj.Image = im;         
end
         

% function cluster = clusterInfluence(cluster)
%         nrC = length(cluster);
%         FullIndex = (1:nrC);
%         for k=1:1:nrC
%             OtherIndex = setdiff(FullIndex,k);
%             %figure;fastrbf_view(warped_cluster(k).uv,warped_cluster(k).Color.Value);
%             for j=1:1:length(OtherIndex)
%                 eval = fastrbf_pointeval(cluster(OtherIndex(j)).ProbRBF,cluster(k).Warped.UV,'messages',0);
%                 %eval = fastrbf_pointeval(warped_cluster(other_index(j)).border.rbf,warped_cluster(k).uv,'messages',0);
%                 KeepIndex = find((cluster(k).Prob-eval.Value)>0);
%                 %keep_index = find((warped_cluster(k).border.Value-eval.Value)>0);
%                 clear eval;
%                 cluster(k).Warped.UV = cluster(k).Warped.UV(:,KeepIndex);
%                 cluster(k).Prob = cluster(k).Prob(KeepIndex);          
%                 cluster(k).Warped.Image = cluster(k).Warped.Image(:,KeepIndex);
%                 clear KeepIndex;
%             end
%         end
% end
         


