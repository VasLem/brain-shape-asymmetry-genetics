function [R] = RotMatrixFromEuler(a,b,c)

    ra=[	1 	 0 	  0 ; 
        0	 cos(a) -sin(a); 
        0 	 sin(a) cos(a)];

    rb=[	cos(b)  0	  sin(b);
        0	  1	  0;
        -sin(b) 0	  cos(b)];
    rc=[	cos(c)  -sin(c)  0;
        sin(c)  cos(c)   0;
        0	  0	     1];

    B=rc*rb*ra;
    R = eye(4,4);
    R(1:3,1:3) = B;
end


