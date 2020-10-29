function out = SampleSizeNPManova(X1,X2,t)
         if nargin < 3, t = 0; end
         nX1 = size(X1,1);
         nX2 = size(X2,1);
         % Getting Subspace centroids & residus
             AvgX1 = mean(X1);
             AvgX2 = mean(X2);
             ResX1 = X1-repmat(AvgX1,nX1,1);
             ResX2 = X2-repmat(AvgX2,nX2,1);
         % Getting within group distances
             PDX1 = sqrt(sum(ResX1.^2,2));
             PDX2 = sqrt(sum(ResX2.^2,2));
         % Getting mean Distances
             out.DispX1 = mean(PDX1);
             out.DispX2 = mean(PDX2);
             X = [PDX1(:);PDX2(:)];
         % generating DistanceMatrix + NP Manova test
           test = oneWayNPManova;
           test.D = squareform(pdist(X,'euclidean'));
           test.n = [nX1 nX2];
           test.t = t;
           perform(test);
           out.Test = test;
end