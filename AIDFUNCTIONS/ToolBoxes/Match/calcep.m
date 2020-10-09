function ep=calcep(a,b)

nrpoints=size(a,2);
ep=0;

for i=1:nrpoints
    tmp=b(:,i)-a(:,i);
    ep=ep + sum(tmp.*tmp);    
end
return;