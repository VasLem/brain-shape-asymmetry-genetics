function ep=calcep_w(a,b,w)

nrpoints=size(a,2);

ep=0;

for i=1:nrpoints
    tmp=w(i)*(b(:,i)-a(:,i));
    ep=ep + sum(tmp.*tmp);    
end

return;