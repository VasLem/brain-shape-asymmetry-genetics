function out=cross(in1,in2)

out(1)=in1(2)*in2(3)-in2(2)*in1(3);
out(2)=in1(3)*in2(1)-in2(3)*in1(1);
out(3)=in1(1)*in2(2)-in2(1)*in1(2);

return;