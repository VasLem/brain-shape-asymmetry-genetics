function [ v ] = plotHierValues3D(values,crit,range)
%F_PLOTPPLOT Circle plot of tree structure, representing p-values
%   Circular plot of values that are structured as a tree. Absolute
%   values are expressed by node size; color represents value signs. 
% 
%   pvalues: 1D input vector where p-values from different levels are
%            concatenated, starting with level 0, etc.
%   Crit: to give special attention to values above the critical value, []
%   when not needed
%   range: if values need to be truncated in a range, do not provide if not
%   wanted, default is auto range

if nargin < 3,
   signvalues = sign(values);
   if sum(signvalues==-1)>0
      range = [-1*max(abs(values)) max(abs(values))];
   else
      range = [min(values) max(values)];
   end
end
if nargin < 2, crit = [];end

if sum(range<0)>0&&sum(range>0)>0
   minrange = 0;
else
   minrange = min(abs(range));
end

%----------- APPREARANCE PARAMETERS -----------%

% Spacing of concentric circles
CIRCLE_SPACING = 12;
Z_SPACING = 10;

% Indicates if concentric circles should be drawn
DRAW_CONCENTRICS = true;

% Min and max size of the markers
MARKER_MIN_SIZE = 1;
MARKER_MAX_SIZE = 10.5;

% Colors of the markers
f = figure;
cm = 'parula';
map = colormap(cm);
close(f);

MARKER_COLOR_POSITIVE = map(end,:);
%MARKER_COLOR_NEGATIVE = map(1,:);

%MARKER_COLOR_POSITIVE = [0 0 0];
%MARKER_COLOR_NEGATIVE = [1 1 1];

% Width of the lines connecting nodes
LINE_WIDTH = 1;

[XS,YS,ZS] = sphere(30);
ZS = -1*ZS;
[FS,VS]  = surf2patch(XS,YS,ZS);
SPHERE = meshObj;
SPHERE.Vertices = VS';
SPHERE.Faces = quad2tri(FS');

% More parameters can be found in code

%-----------------------------------------------%


%values = -log10(values);
%crit = -log10(crit);

values(values<range(1)) = range(1);
values(values>range(2)) = range(2);

% Size of marker expresses absolute p-values

marker_sizes = abs(values)-minrange;
marker_sizes = (marker_sizes/max(abs(range)))*(MARKER_MAX_SIZE-MARKER_MIN_SIZE)+MARKER_MIN_SIZE;

%marker_sizes = abs(values)-min(abs(values));
%marker_sizes = (marker_sizes/max(marker_sizes))*(MARKER_MAX_SIZE-MARKER_MIN_SIZE)+MARKER_MIN_SIZE;


v = viewer3DObj; v.BackgroundColor = [1 1 1];v.AxesVisible = false;v.AxesWallColor = [1 1 1];
v.AxesXColor = [0 0 0];
v.AxesYColor = [0 0 0];
v.AxesZColor = [0 0 0];

% v.SceneLightVisible = true;
% v.SceneLightLinked = true;

hold(v.RenderAxes,'on');
nL = (log2(length(values)+1)-1);
previous_depth = nL*Z_SPACING;
circle_centers  = [0;0;previous_depth];
previous_points = [0;0;previous_depth];
for L = 1:(log2(length(values)+1)-1)
    %L = 1;
    point_count     = 2^L;
    point_list      = [];
    angle_increment = 360/point_count;
    
    % Calculate points at current level
    for P = 1:point_count
        a = angle_increment*(P-1) + angle_increment/2;
        a = a+90;
        point_list(1:3,P) = [cosd(a) -sind(a) 0; sind(a) cosd(a) 0; 0 0 1]*[CIRCLE_SPACING*L;0;1]; %#ok<*AGROW>
    end
    point_list(3,:) = previous_depth-Z_SPACING;
    
    % Draw concentric if needed
    if DRAW_CONCENTRICS
        ang=0:0.01:2*pi; 
        xp=L*CIRCLE_SPACING*cos(ang);
        yp=L*CIRCLE_SPACING*sin(ang);
        zp = (previous_depth-Z_SPACING)*ones(1,length(ang));
        plot3(v.RenderAxes,xp,yp,zp,'Color','k','LineStyle',':','LineWidth',LINE_WIDTH);
    end
    
    % Draw connections from previous to current
    for P = 1:2^(L-1)
        plot3(v.RenderAxes,[previous_points(1,P),point_list(1,(P-1)*2+1)],....
                           [previous_points(2,P),point_list(2,(P-1)*2+1)],...
                           [previous_points(3,P),point_list(3,(P-1)*2+1)],'k','linewidth',LINE_WIDTH);
        plot3(v.RenderAxes,[previous_points(1,P),point_list(1,(P-1)*2+2)],...
                           [previous_points(2,P),point_list(2,(P-1)*2+2)],...
                           [previous_points(3,P),point_list(3,(P-1)*2+2)],'k','linewidth',LINE_WIDTH);
    end 
    % Current points become previous; centers are stored
    previous_points = point_list;
    previous_depth = previous_depth-Z_SPACING;
    circle_centers(:,end+1:end+size(point_list,2)) = point_list(1:3,:);
end

% Plot nodes on top
for n=1:size(circle_centers,2)
    %n = 1;
    scan = clone(SPHERE);
    scan.Vertices = (marker_sizes(n)/2*scan.Vertices)+repmat(circle_centers(:,n),1,scan.nrV);
    scan.Axes = v.RenderAxes;
    
    ang=0:0.01:2*pi; 
    xp=marker_sizes(n)/2*cos(ang)+circle_centers(1,n);
    yp=marker_sizes(n)/2*sin(ang)+circle_centers(2,n);
    zp = circle_centers(3,n)*ones(1,length(ang));
    if ~isempty(crit)
        if values(n) > crit
            scan.Material = 'Dull';
            %scan.Material = 'Default';
            plot3(v.RenderAxes,xp,yp,zp,'Color',MARKER_COLOR_POSITIVE,'LineStyle','-','LineWidth',LINE_WIDTH*5);
        else
            scan.Material = 'Dull';
            %scan.Material = 'Default';
            %plot3(v.RenderAxes,xp,yp,zp,'Color',MARKER_COLOR_NEGATIVE,'LineStyle','-','LineWidth',LINE_WIDTH*3);
        end
    else
        scan.Material = 'Dull';
    end
    scan.Value = values(n)*ones(1,scan.nrV);
    scan.ColorMode = 'Indexed';
    scan.Visible = true;
    scan.Selected = true;
end
colormap(v.RenderAxes,cm);colorbar('peer',v.RenderAxes);
set(v.RenderAxes,'clim',range);
%v.CamProjection = 'perspective';
% Set the axis appropriately
max_space = MARKER_MAX_SIZE/2;
xlim([min(point_list(1,:))-max_space,max(point_list(1,:))+max_space]);
ylim([min(point_list(2,:))-max_space,max(point_list(2,:))+max_space]);
hold off; axis equal;drawnow;
end

