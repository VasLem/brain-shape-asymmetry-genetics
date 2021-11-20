function res = assignPvalueSignificance(pValues,smallestPvalue)
   if nargin == 1, smallestPvalue=0.0001; end
   thresholds = 10.^-(log10(1/smallestPvalue)-(2:-1:0));
    thresholds = [thresholds * 5; thresholds];
    thresholds = thresholds(:)';
   res = zeros(size(pValues));
   c = 1;
   for t=thresholds
    res(pValues <= t) = c;
    c = c + 1;
   end
end