function [R]=corr2mat(X,Y)
%USAGE [R]=corr2mat(X,Y)
% Correlation between 2 conformable matrices
% X is I*J, Y is I*K
% if # rows is different -> error
% R is J*K 
% Herve Abdi February 2006
if nargin==1;
    R=corrcoef(X);
    return
end
[nix,njx]=size(X);
[niy,njy]=size(Y);
if nix~=niy;
    error(['corr2mat -> Input matrices non conformable']);
end
R=corrcoef([X Y]);
R=R(1:njx,njx+1:njx+njy);

