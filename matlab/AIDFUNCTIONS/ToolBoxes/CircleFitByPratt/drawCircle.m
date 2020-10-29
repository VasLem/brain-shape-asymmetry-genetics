function drawCircle(par,nseg)

theta = 0 : (2 * pi / nseg) : (2 * pi);
pline_x = par(3) * cos(theta) + par(1);
pline_y = par(3) * sin(theta) + par(2);

plot(pline_x,pline_y);

end
