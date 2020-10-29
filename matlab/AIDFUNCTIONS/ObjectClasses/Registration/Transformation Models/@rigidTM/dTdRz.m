function out = dTdRz(obj)
        Rx = obj.P(1);Ry = obj.P(2);Rz = obj.P(3);
        out = zeros(3,3);
        out(1,1) = -sin(Rz)*cos(Ry); 
        out(1,2) = -cos(Rz)*cos(Rx)-sin(Rz)*sin(Ry)*sin(Rx);
        out(1,3) = cos(Rz)*sin(Rx)-sin(Rz)*sin(Ry)*cos(Rx);
        out(2,1) = cos(Rz)*cos(Ry);
        out(2,2) = -sin(Rz)*cos(Rx)+cos(Rz)*sin(Ry)*sin(Rx);
        out(2,3) = sin(Rz)*sin(Rx)+cos(Rz)*sin(Ry)*cos(Rx);
end
%out(4,4) = 1;