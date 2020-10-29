function out = dTdRy(obj)
        Rx = obj.P(1);Ry = obj.P(2);Rz = obj.P(3);
        out = zeros(3,3);
        out(1,1) = -cos(Rz)*sin(Ry);
        out(1,2) = cos(Rz)*cos(Ry)*sin(Rx);
        out(1,3) = cos(Rz)*cos(Ry)*cos(Rx);
        out(2,1) = -sin(Rz)*sin(Ry);
        out(2,2) = sin(Rz)*cos(Ry)*sin(Rx);
        out(2,3) = sin(Rz)*cos(Ry)*cos(Rx);
        out(3,1) = -cos(Ry);
        out(3,2) = -sin(Ry)*sin(Rx);
        out(3,3) = -sin(Ry)*cos(Rx);
end
