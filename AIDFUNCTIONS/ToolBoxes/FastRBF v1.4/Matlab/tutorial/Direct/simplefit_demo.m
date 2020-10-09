% Demonstration file for a simple RBF biharmonic spline 
% fit in 3D
%
%
n=100;
cent=rand(100,3);
f=rand(100,1);

[A]=spline_direct(cent);

coeff=spline_fitit(cent,f);
fprintf('Value to be interpolated  %15.10e \n',f(1));
fprintf('Value of RBF interpolant  %15.10e \n',...
       spline_eval(cent,coeff,cent(1,:)));
