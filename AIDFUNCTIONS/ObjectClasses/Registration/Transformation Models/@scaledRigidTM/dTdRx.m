function out = dTdRx(obj)
    Rx = obj.P(1);Ry = obj.P(2);Rz = obj.P(3);
    out = zeros(3,3);
    out(1,2) = (sin(Rz)*sin(Rx)+cos(Rz)*sin(Ry)*cos(Rx));
    out(1,3) = (sin(Rz)*cos(Rx)-cos(Rz)*sin(Ry)*sin(Rx));
    out(2,2) = (-cos(Rz)*sin(Rx)+sin(Rz)*sin(Ry)*cos(Rx));
    out(2,3) = (-cos(Rz)*cos(Rx)-sin(Rz)*sin(Ry)*sin(Rx));
    out(3,2) = cos(Ry)*cos(Rx);
    out(3,3) = -cos(Ry)*sin(Rx);
end