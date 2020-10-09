function [ output_args ] = f_plotPValues( pvalues, pCrit )
%F_PLOTPPLOT Circle plot of tree structure, representing p-values
%   Circular plot of p-values that are structured as a tree. Absolute
%   p-values are expressed by node size; color represents p-value signs. 
% 
%   pvalues: 1D input vector where p-values from different levels are
%            concatenated, starting with level 0, etc.


%----------- APPREARANCE PARAMETERS -----------%

% Spacing of concentric circles
CIRCLE_SPACING = 12;

% Indicates if concentric circles should be drawn
DRAW_CONCENTRICS = true;

% Min and max size of the markers
MARKER_MIN_SIZE = 3;
MARKER_MAX_SIZE = 12;

% Colors of the markers
MARKER_COLOR_POSITIVE = 'b';
MARKER_COLOR_NEGATIVE = 'r';

% Width of the lines connecting nodes
LINE_WIDTH = 2;

% More parameters can be found in code

%-----------------------------------------------%

pvalues = -log10(pvalues);
if nargin < 2, pCrit = 1;end
pCrit = -log10(pCrit);
% Size of marker expresses absolute p-values
marker_sizes = abs(pvalues)-min(abs(pvalues));
marker_sizes = (marker_sizes/max(marker_sizes))*(MARKER_MAX_SIZE-MARKER_MIN_SIZE)+MARKER_MIN_SIZE;

figure; hold on;
circle_centers  = [0;0];
previous_points = [0;0];
marker_index    = 1;
for L = 1:(log2(length(pvalues)+1)-1)
    point_count     = 2^L;
    point_list      = [];
    angle_increment = 360/point_count;
    
    % Calculate points at current level
    for P = 1:point_count
        a = angle_increment*(P-1) + angle_increment/2;
        point_list(1:3,P) = [cosd(a) -sind(a) 0; sind(a) cosd(a) 0; 0 0 1]*[CIRCLE_SPACING*L;0;1];
    end
    
    % Draw concentric if needed
    if DRAW_CONCENTRICS
        viscircles([0,0],L*CIRCLE_SPACING,'Color','k','LineStyle',':','LineWidth',2);
    end
    
    % Draw connections from previous to current
    for P = 1:2^(L-1)
        plot([previous_points(1,P),point_list(1,(P-1)*2+1)],[previous_points(2,P),point_list(2,(P-1)*2+1)],'k','linewidth',LINE_WIDTH);
        plot([previous_points(1,P),point_list(1,(P-1)*2+2)],[previous_points(2,P),point_list(2,(P-1)*2+2)],'k','linewidth',LINE_WIDTH);
        marker_index = marker_index + 2;
    end
    
    % Current points become previous; centers are stored
    previous_points = point_list;
    circle_centers(:,end+1:end+size(point_list,2)) = point_list(1:2,:);
end

% Plot nodes on top
for n=1:size(circle_centers,2)
    x = circle_centers(1,n)-marker_sizes(n)/2;
    y = circle_centers(2,n)-marker_sizes(n)/2;
    if pvalues(n) > pCrit
        rectangle('Position',[x,y,marker_sizes(n),marker_sizes(n)],'Curvature',[1 1],'FaceColor',MARKER_COLOR_POSITIVE);
    else
        rectangle('Position',[x,y,marker_sizes(n),marker_sizes(n)],'Curvature',[1 1],'FaceColor',MARKER_COLOR_NEGATIVE);
    end
end

% Set the axis appropriately
max_space = MARKER_MAX_SIZE/2;
xlim([min(point_list(1,:))-max_space,max(point_list(1,:))+max_space]);
ylim([min(point_list(2,:))-max_space,max(point_list(2,:))+max_space]);
hold off; axis equal;
end

