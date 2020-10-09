function z = P2Z(p)
    %z1 = -sqrt(2) * erfcinv((1-p/2)*2);
    z = abs(-sqrt(2) * erfcinv((p/2)*2));
end