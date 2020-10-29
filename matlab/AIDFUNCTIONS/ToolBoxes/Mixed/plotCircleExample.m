close all; clc;

CIRCLE_SPACING = 10;
CIRCLE_LEVELS  = 5;
MARKER_SIZE    = 60;
LINE_WIDTH     = 3;
ADD_CIRCLES    = true;

figure; hold on;
previous_points = [0;0];
for C = 1:CIRCLE_LEVELS
    if ADD_CIRCLES
        viscircles([0,0],C*CIRCLE_SPACING);
    end
    point_count = 2^C;
    point_list  = [];
    start_point = [CIRCLE_SPACING*C;0;1];
    angle_increment = 360/point_count;
    for P = 1:point_count
        a = angle_increment*(P-1) + angle_increment/2;
        point_list(1:3,P) = [cosd(a) -sind(a) 0; sind(a) cosd(a) 0; 0 0 1]*start_point;
    end
    for P = 1:2^(C-1)
        plot([previous_points(1,P),point_list(1,(P-1)*2+1)],[previous_points(2,P),point_list(2,(P-1)*2+1)],'o-r','linewidth',LINE_WIDTH,'markersize',MARKER_SIZE,'markerfacecolor','g');
        plot([previous_points(1,P),point_list(1,(P-1)*2+2)],[previous_points(2,P),point_list(2,(P-1)*2+2)],'o-r','linewidth',LINE_WIDTH,'markersize',MARKER_SIZE,'markerfacecolor','g');
    end
    previous_points = point_list;
end
hold off; axis equal;

load EuropeanHierarchical_20112015;
load RefScan;
RefScan.Material = 'Dull';
RefScan.Value = ones(1,RefScan.nrV);
RefScan.ColorMode = 'Indexed';
scan = clone(RefScan);
% FOR EACH LEVEL, ALL CLUSTERS
v = viewer3DObj;
pause;
for l=1:1:size(label,2)
   % l=2;
   for i=1:2^(l-1)
       % i=1
      scan = clone(RefScan);
      scan.Value = double((label(l,:)==i));
      pos = 10*rand(3,1);
      scan.Vertices = scan.Vertices-repmat(pos,1,scan.nrV);
      scan.Axes = v.RenderAxes;
      scan.Visible = true;
      pause;
   end
end

