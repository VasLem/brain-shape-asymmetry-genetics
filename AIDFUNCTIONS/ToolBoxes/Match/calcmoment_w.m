function M=calcmoment_w(a,b,w)

nrpoints=size(a,2);
M=zeros(1,3);
for i=1:nrpoints
    M=M+w(i)*cross2(b(:,i),a(:,i)-b(:,i));
end
return;