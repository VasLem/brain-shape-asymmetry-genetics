function plotline(F,G,ax1,ax2,mark)
% In development Jan 99
%USAGE: plotligne(F,G,ax1,ax2,mark)
% plot point with lines between them 
% F coor of mat 1
% G coor of mat 2
% ax1, ax2 the axis from F and G
% mark: the definition of the line
% eg '-r' '-.m'
[nI_F,nK_F]=size(F);
if exist('ax1')==0;ax1=1;end
if exist('ax2')==0;ax2=2;end
if exist('mark')==0;mark='-';end

for i=1:nI_F;
   x=[F(i,ax1),G(i,ax1)];
   y=[F(i,ax2),G(i,ax2)];
   plot(x,y,mark)
end
