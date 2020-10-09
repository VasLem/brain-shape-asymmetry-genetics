function [A]=direct(cent)
% Form the (n+dim+1)*(n+dim+1) matrix corresponding to 
% RBF interpolation with basic function phi and linears.
%
% The centres are assumed given in the n by dim array cent.
% phi is assumed given as a function of r. It is coded in 
% the Matlab function phi. m
%
% Syntax  [A]=direct(cent)
%
% Input          
%         cent  n*dim array  coordinates of centers
%               of reals     
% Output  A     (n+dim+1)*   Symmetric matrix of the linear
%               (n+dim+1)    system obtained if we solve
%               array of     for the radial basis function
%                            interpolant directly.
%
% Write the matrix A in the form
%             B    P
%     A   =  
%             P^t  O
% where P is the polynomial bit.
%
[n dim]=size(cent);
%A=zeros(n,n);
parfor i=1:n
    centfor = cent;
    Afor = zeros(1,i);
    for j=1:i
        r=norm(cent(i,:)-centfor(j,:));
        temp=phi(r);
        Afor(1,j)=temp;
        %A(j,i)=temp;
    end
    A{i} = Afor;
end
AN = zeros(n,n);
for i=1:1:n
    AN(i,1:i) = AN{i};
end
keyboard;
%
% Now the polynomial part
%
P=[ones(n,1) cent];
A = [ A      P
      P' zeros(dim+1,dim+1)];

