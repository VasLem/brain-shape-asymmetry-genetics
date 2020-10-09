function M=calcmoment(a,b)

nrpoints=size(a,2);
M=zeros(1,3);
for i=1:nrpoints
    M=M+cross2(b(:,i),a(:,i)-b(:,i));
end
return;