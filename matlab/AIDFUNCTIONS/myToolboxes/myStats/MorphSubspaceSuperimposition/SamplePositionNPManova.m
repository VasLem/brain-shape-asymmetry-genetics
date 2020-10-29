function out = SamplePositionNPManova(X1,X2,t)
         if nargin < 3, t = 0; end
         nX1 = size(X1,1);
         nX2 = size(X2,1);
         % Getting Centroids
           AvgX1 = mean(X1);
           AvgX2 = mean(X2);
         % generating Effect output  
           out.AvgX1 = AvgX1;
           out.AvgX2 = AvgX2;
           out.Difference = AvgX1-AvgX2;
           out.Distance = sqrt(sum(out.Difference.^2));
         % generating DistanceMatrix + NP Manova test
           test = oneWayNPManova;
           test.D = squareform(pdist([X1;X2],'euclidean'));
           test.n = [nX1 nX2];
           test.t = t;
           perform(test);
           out.Test = test;
end