function res = assignPvalueSignificance(pValues,thresholds)
   if nargin == 1, thresholds = [0.05, 0.01, 0.005, 0.0001]; end
   res = zeros(size(pValues));
   c = 1;
   for t=thresholds
    res(pValues <= t) = c;
    c = c + 1;
   end
end