function [coeff]=fitit(cent,f)
% Function to find the radial basis function s consisting
% of a linear plus a sum of shifts of \phi( | x | )
% interpolating to data
% f_i at the point cent(i,:) for 1 \leq i \leq n.
%
% Syntax [coeff]=fitit(cent,f)
%
% Input
%
%      cent    n by dim array of centres
%      f       n by 1 vector of values at centres
%
% Output
%      coeff   (n+dim+1) vector of coefficients with the coefficients
%              of 1, x and y (2D) or 1, x, y and z (3D) last.
%

A=direct(cent);
[n dim] =size(cent);
f=f(:); % Make sure its a column
% add the zeros for the polynomial part
% at the end of the column
f=[f ;zeros(dim+1,1)]; 
coeff=A\f;
