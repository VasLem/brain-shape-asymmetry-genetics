function ret = checkLSize(X1, X2, lSizes)
assert(max(lSizes)<= size(X1,2)/3);
numExp = length(lSizes);
for i=1:numExp
    inputSize = lSizes(i);
    subsample_vec = 3 * (1: ceil((size(X1,2)/(3 * inputSize))):(size(X1,2)/3));
    subsample_vec = ([ (subsample_vec-2)' (subsample_vec -1)' subsample_vec'])';
    subsample_vec = subsample_vec(:);
    ret(i)  = ProcrustesAnova2WayAsymmetryMEM(X1(:, subsample_vec,:),X2(:,subsample_vec,:),200);
end