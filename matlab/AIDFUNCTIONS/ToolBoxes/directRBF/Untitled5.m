% Demonstration file for a simple RBF biharmonic spline 
% fit in 3D
%
%
n=10000;
cent=rand(n,3);
f=rand(n,1);


tic
coeff=fitit(cent,f);
toc;




fprintf('Value to be interpolated  %15.10e \n',f(1));
fprintf('Value of RBF interpolant  %15.10e \n',...
       eval_direct(cent,coeff,cent(1,:)));

