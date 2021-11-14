function ret = checkRepsNum(X1, X2, reps,factor)
assert(round(size(X1,3)/max(reps)) ==size(X1,3)/max(reps))
numExp = length(reps);

for i=1:numExp
    inputReps = reps(i);
    subsample_vec = 1: ceil((size(X1,3)/inputReps)):size(X1,3);
    ret(i) = ProcrustesAnova2WayAsymmetryMEM(X1(:,:,subsample_vec),X2(:,:,subsample_vec),1000,factor);
end
end