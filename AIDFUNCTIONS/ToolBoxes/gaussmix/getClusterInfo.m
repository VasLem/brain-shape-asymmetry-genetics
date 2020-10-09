function [cluster] = getClusterInfo(mtrs,points)
    nr_cluster = size(mtrs.cluster,2);
    for k=1:1:nr_cluster
            cluster(k).mu = mtrs.cluster(k).mu;
            cluster(k).cov = mtrs.cluster(k).R;
            cluster(k).max_prob = mvnpdf(cluster(k).mu',cluster(k).mu',cluster(k).cov);
            y  = mvnpdf(points',cluster(k).mu',cluster(k).cov);
            cluster(k).y = y/cluster(k).max_prob;        
            cluster(k).min_prob = min(y);
    end    
end