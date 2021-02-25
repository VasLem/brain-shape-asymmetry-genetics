function ret = checkSize(X1, X2, sizes)
assert(max(sizes)<= size(X1,1))
numExp = length(sizes);
ret = [];

parfor i=1:numExp
    inputSize = sizes(i);
    subsample_vec = 1: ceil((size(X1,1)/inputSize)):size(X1,1);
    append(ret, ProcrustesAnova2WayAsymmetryMEM(X1(subsample_vec,:,:),X2(subsample_vec,:,:),100));
end
end