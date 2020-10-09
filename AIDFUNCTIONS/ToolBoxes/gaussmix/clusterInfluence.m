function cluster = clusterInfluence(cluster)
        nrC = length(cluster);
        FullIndex = (1:nrC);
        for k=1:1:nrC
            OtherIndex = setdiff(FullIndex,k);
            %figure;fastrbf_view(warped_cluster(k).uv,warped_cluster(k).Color.Value);
            for j=1:1:length(OtherIndex)
                eval = fastrbf_pointeval(cluster(OtherIndex(j)).ProbRBF,cluster(k).Warped.UV,'messages',0);
                %eval = fastrbf_pointeval(warped_cluster(other_index(j)).border.rbf,warped_cluster(k).uv,'messages',0);
                KeepIndex = find((cluster(k).Prob-eval.Value)>-0.01);
                %keep_index = find((warped_cluster(k).border.Value-eval.Value)>0);
                clear eval;
                cluster(k).Warped.UV = cluster(k).Warped.UV(:,KeepIndex);
                cluster(k).Prob = cluster(k).Prob(KeepIndex);          
                cluster(k).Warped.Image = cluster(k).Warped.Image(:,KeepIndex);
                clear KeepIndex;
            end
        end
end