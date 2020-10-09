function [out,info] = NormalizeExpressionSeverity(ass,regionindex,sev)
         regdistances = ass.DistanceMap(regionindex);
         origsev = sqrt(mean((regdistances.^2)));
         [xout,opterror] = fminsearch(@(X) errorfunction(X,ass,regionindex,sev),0);
         tmp = clone(ass.Scan);
         tmp.Vertices = xout*ass.Scan.Vertices + (1-xout)*ass.Norm.Vertices;
         out = assessment;
         out.Scan = clone(ass.Scan);
         out.Norm = tmp;
         out.DistanceRange = ass.DistanceRange;
         update(out);
         info.X = xout;
         info.opterror = opterror;
         info.origsev = origsev;
         
end

function err = errorfunction(X,ass,regionindex,sev)
         tmp = clone(ass.Scan);
         tmp.Vertices = X*ass.Scan.Vertices + (1-X)*ass.Norm.Vertices;
         dist = vDistances(ass.Scan,tmp);
         err = abs(sev-sqrt(mean((dist(regionindex).^2))));
end