function [ v ] = plotHierFaces3D(RefScan,label)
%F_PLOTPPLOT Circle plot of tree structure, representing p-values
%   Circular plot of p-values that are structured as a tree. Absolute
%   p-values are expressed by node size; color represents p-value signs. 
% 
%   pvalues: 1D input vector where p-values from different levels are
%            concatenated, starting with level 0, etc.


%----------- APPREARANCE PARAMETERS -----------%

% Spacing of concentric circles
CIRCLE_SPACING = 12;
Z_SPACING = 10;

% Indicates if concentric circles should be drawn
DRAW_CONCENTRICS = true;

% Min and max size of the markers
MARKER_MIN_SIZE = 3;
MARKER_MAX_SIZE = 15;

% Colors of the markers
MARKER_COLOR_POSITIVE = 'b';
MARKER_COLOR_NEGATIVE = 'r';

% Width of the lines connecting nodes
LINE_WIDTH = 1;

[XS,YS,ZS] = sphere(40);
ZS = -1*ZS;
[FS,VS]  = surf2patch(XS,YS,ZS);
SPHERE = meshObj;
SPHERE.Vertices = VS';
SPHERE.Faces = quad2tri(FS');

RefScan = clone(RefScan);
theta = -pi/4;
ROT = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
RefScan.Vertices = ROT*RefScan.Vertices;
% More parameters can be found in code

%-----------------------------------------------%


BK = [1 1 1];
FG = [0 0 0];


v = viewer3DObj; v.BackgroundColor = BK;v.AxesVisible = false;

v.SceneLightVisible = true;
v.SceneLightLinked = true;


hold on;
%nL = (log2(length(values)+1)-1);
nL = size(label,1)-1;
previous_depth = nL*Z_SPACING;
circle_centers  = [0;0;previous_depth];
previous_points = [0;0;previous_depth];
for L = 1:nL
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
        plot3(v.RenderAxes,xp,yp,zp,'Color',FG,'LineStyle',':','LineWidth',LINE_WIDTH);
    end
    
    % Draw connections from previous to current
    for P = 1:2^(L-1)
        plot3(v.RenderAxes,[previous_points(1,P),point_list(1,(P-1)*2+1)],....
                           [previous_points(2,P),point_list(2,(P-1)*2+1)],...
                           [previous_points(3,P),point_list(3,(P-1)*2+1)],'Color',FG,'linewidth',LINE_WIDTH);
        plot3(v.RenderAxes,[previous_points(1,P),point_list(1,(P-1)*2+2)],...
                           [previous_points(2,P),point_list(2,(P-1)*2+2)],...
                           [previous_points(3,P),point_list(3,(P-1)*2+2)],'Color',FG,'linewidth',LINE_WIDTH);
    end 
    % Current points become previous; centers are stored
    previous_points = point_list;
    previous_depth = previous_depth-Z_SPACING;
    circle_centers(:,end+1:end+size(point_list,2)) = point_list(1:3,:);
end

levels = getLevels(nL);
clusters = getClusters(nL);

% Plot nodes on top
for n=1:size(circle_centers,2)
    %n = 1;
    scan = clone(RefScan);
    scan.Vertices = 4*scan.Vertices+repmat(circle_centers(:,n),1,scan.nrV);
    scan.Material = 'Dull';
    tmplev = levels(n);
    tmpcl = clusters(n);
    scan.Value = double((label(tmplev,:)==tmpcl));
    scan.Axes = v.RenderAxes;
    scan.ColorMode = 'Indexed';
    scan.Visible = true;
    scan.Selected = true;
end
colormap(v.RenderAxes,'summer');%colorbar('peer',v.RenderAxes);
v.CamProjection = 'perspective';
% Set the axis appropriately
max_space = MARKER_MAX_SIZE/2;
xlim([min(point_list(1,:))-max_space,max(point_list(1,:))+max_space]);
ylim([min(point_list(2,:))-max_space,max(point_list(2,:))+max_space]);
hold off; axis equal;
end


function out = getLevels(L)
            out = [];
            for i=1:1:L+1
               out = [out i*ones(1,2^(i-1))];  %#ok<AGROW>
            end
end
function out = getClusters(L)
           out = [];
           for i=1:1:L+1
               out = [out 1:2^(i-1)]; %#ok<AGROW>
           end
end

