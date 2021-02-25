function ret = checkSize(X1, X2, sizes,cntFromBeg)
if nargin < 4
    cntFromBeg=0;
end
assert(max(sizes)<= size(X1,1))
numExp = length(sizes);


for i=1:numExp
    inputSize = sizes(i);
    if ~cntFromBeg
        
        subsample_vec = 1: ceil((size(X1,1)/inputSize)):size(X1,1);
    else
        subsample_vec = 1:inputSize;
    end
    ret(i) = ProcrustesAnova2WayAsymmetryMEM(X1(subsample_vec,:,:),X2(subsample_vec,:,:),100);
end
end

