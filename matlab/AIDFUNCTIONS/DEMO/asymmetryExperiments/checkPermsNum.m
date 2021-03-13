function ret = checkPermsNum(X1, X2, perms)
numExp = length(perms);
for i=1:numExp
    p = perms(i);
    ret(i)  = ProcrustesAnova2WayAsymmetryMEM(X1,X2,p);
end
end

