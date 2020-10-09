function plotxyrange(F,axe1,axe2,titre,noms,range);
% % ***** This is a Test  Version *******
% January  1999 Herve Abdi
% Usage plotxyabs(F,axe1,axe2,title,names,range);
% plotxy plots a MDS or PCA or CA graph of component #'Axe1' vs #'Axe2'
% F is a matrix of coordinates
% axe1 is the Horizontal Axis
% axe2 is the Vertical Axis
% title will be the title of the graph
% range is a 1*4 vector giving
% the [minx,maxx,miny,maxy] values (def as plotxyha)
% Axes are labelled 'Principal Component number'
% names give the names of the points to plot (def=numbers)
% see also plotxyha , plotxyha2, plotxyabs
% difference in dealing with the names and the first axis
% plotxyabs plot the First axis as starting at 0.


[nrow,ncol]=size(F);
if exist('noms')==0;
   noms{nrow,1}=[];
   for k=1:nrow;noms{k,1}=int2str(k);end
end
if isa(noms,'char')==1;
   nomdenom=noms;clear noms;
   noms{nrow,1}=[];
   for k=1:nrow;noms{k,1}=nomdenom(k,:);end
end

if exist('range')==0;
 minx=min(F(:,axe1));
 maxx=max(F(:,axe1));
 miny=min(F(:,axe2));
 maxy=max(F(:,axe2));
range=[minx,maxx,miny,maxy]; 
%[minx,maxx,miny,maxy]=range(:);
else
minx=range(1);maxx=range(2);
miny=range(3);maxy=range(4);
end
hold off; clf;hold on;
 rangex=maxx-minx;epx=rangex/10;
 rangey=maxy-miny;epy=rangey/10; axis('equal');
 axis([minx-epx,maxx+epx,miny-epy,maxy+epy]) ;
 %axis('equal');
%axis;
plot ( F(:,axe1),F(:,axe2 ),'.');
label=' Principal Component Number ';
labelx=[label, num2str(axe1) ];
labely=[label, num2str(axe2) ];
xlabel (labelx);
ylabel (labely );
plot([minx-epx,maxx+epx],[0,0] ,'b');
% hold
plot ([0,0],[miny-epy,maxy+epy],'b');
for i=1:nrow,
  text(F(i,axe1),F(i,axe2),noms{i,1});
end;
title(titre);

