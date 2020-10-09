function [ v,circle_centers ] = plotPaperCombined(values,varargin)
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

%------------- SCREENING INPUT ----------------%
    Input = find(strcmpi(varargin, 'crit'));
    if isempty(Input)
       crit = [];
    else
       crit = varargin{Input+1};
    end
    Input = find(strcmpi(varargin, 'range'));
    if isempty(Input)
       signvalues = sign(values);
       if sum(signvalues==-1)>0
          range = [-1*max(abs(values)) max(abs(values))];
       else
          range = [min(values) max(values)];
       end
    else
       range = varargin{Input+1};
    end
    Input = find(strcmpi(varargin, 'light'));
    if isempty(Input)
       light = false;
    else
       light = varargin{Input+1};
    end
    Input = find(strcmpi(varargin, 'peer'));
    if isempty(Input)
       v = viewer3DObj; v.BackgroundColor = [1 1 1];v.AxesVisible = false;v.AxesWallColor = [1 1 1];
       v.AxesXColor = [0 0 0];v.AxesYColor = [0 0 0];v.AxesZColor = [0 0 0];
       plotax = v.RenderAxes;
    else
       v = varargin{Input+1};
       switch v.Type
           case 'axes'
               plotax = v;
               hold(plotax,'on');
           case 'viewer3DObj'
               plotax = v.RenderAxes;
           case 'figure'
               plotax = v.CurrentAxes;
               hold(plotax,'off');
       end
    end
    %cla(plotax);
    Input = find(strcmpi(varargin, 'labels'));
    if isempty(Input)
       labels = {};
    else
       labels = varargin{Input+1};
       if isnumeric(labels)
          tmp = cell(1,length(labels));
          for i=1:length(tmp)
              tmp{i} = num2str(labels(i));
          end
          labels = tmp; 
       end
    end
    Input = find(strcmpi(varargin, 'title'));
    if isempty(Input)
       titlestr = [];
    else
       titlestr = varargin{Input+1};
    end
    Input = find(strcmpi(varargin, 'circlespacing'));
    if isempty(Input)
       CIRCLE_SPACING = 12;
    else
       CIRCLE_SPACING = varargin{Input+1};
    end
    Input = find(strcmpi(varargin, 'zspacing'));
    if isempty(Input)
       Z_SPACING = 10;
    else
       Z_SPACING = varargin{Input+1};
    end
    Input = find(strcmpi(varargin, 'concentrics'));
    if isempty(Input)
       DRAW_CONCENTRICS = true;
    else
       DRAW_CONCENTRICS = varargin{Input+1};
    end
    Input = find(strcmpi(varargin, 'minsize'));
    if isempty(Input)
       MARKER_MIN_SIZE = 1;
    else
       MARKER_MIN_SIZE = varargin{Input+1};
    end
    Input = find(strcmpi(varargin, 'maxsize'));
    if isempty(Input)
       MARKER_MAX_SIZE = 10.5;
    else
       MARKER_MAX_SIZE = varargin{Input+1};
    end
    Input = find(strcmpi(varargin, 'linewidth'));
    if isempty(Input)
       LINE_WIDTH = 1;
    else
       LINE_WIDTH = varargin{Input+1};
    end
    Input = find(strcmpi(varargin, 'colormap'));
    if isempty(Input)
       cm = 'parula';
    else
       cm = varargin{Input+1};
    end
    %----------- APPREARANCE PARAMETERS -----------%

    % Colors of the markers
    f = figure;map = colormap(cm);close(f);
    MARKER_COLOR_POSITIVE = map(end,:);

    % CREATING UNIT SPHERE
    [XS,YS,ZS] = sphere(30);
    ZS = -1*ZS;
    [FS,VS]  = surf2patch(XS,YS,ZS);
    SPHERE = meshObj;
    SPHERE.Vertices = VS';
    SPHERE.Faces = quad2tri(FS');

    % More parameters can be found in code

    %-----------------------------------------------%

    if sum(range<0)>0&&sum(range>0)>0
       minrange = 0;
    else
       minrange = min(abs(range));
    end
    values(values<range(1)) = range(1);
    values(values>range(2)) = range(2);
    % Size of marker expresses absolute values
    marker_sizes = abs(values)-minrange;
    marker_sizes = (marker_sizes/max(abs(range)))*(MARKER_MAX_SIZE-MARKER_MIN_SIZE)+MARKER_MIN_SIZE;

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
            plot3(plotax,xp,yp,zp,'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',LINE_WIDTH);
            hold(plotax,'on');
        end

        % Draw connections from previous to current
        for P = 1:2^(L-1)
            plot3(plotax,[previous_points(1,P),point_list(1,(P-1)*2+1)],....
                               [previous_points(2,P),point_list(2,(P-1)*2+1)],...
                               [previous_points(3,P),point_list(3,(P-1)*2+1)],'Color',[0.7 0.7 0.7],'linewidth',LINE_WIDTH);
            plot3(plotax,[previous_points(1,P),point_list(1,(P-1)*2+2)],...
                               [previous_points(2,P),point_list(2,(P-1)*2+2)],...
                               [previous_points(3,P),point_list(3,(P-1)*2+2)],'Color',[0.7 0.7 0.7],'linewidth',LINE_WIDTH);
        end 
        % Current points become previous; centers are stored
        previous_points = point_list;
        previous_depth = previous_depth-Z_SPACING;
        circle_centers(:,end+1:end+size(point_list,2)) = point_list(1:3,:);
    end

    
    
    return;    
    
    
    
    
    % Plot nodes on top
    for n=1:size(circle_centers,2)
        %n = 1;
        scan = clone(SPHERE);
        scan.Vertices = (marker_sizes(n)/2*scan.Vertices)+repmat(circle_centers(:,n),1,scan.nrV);
        scan.Axes = plotax;

        ang=0:0.01:2*pi; 
        xp=marker_sizes(n)/2*cos(ang)+circle_centers(1,n);
        yp=marker_sizes(n)/2*sin(ang)+circle_centers(2,n);
        zp = circle_centers(3,n)*ones(1,length(ang));
        if ~isempty(crit)
            if values(n) > crit
                scan.Material = 'Dull';
                plot3(plotax,xp,yp,zp,'Color',MARKER_COLOR_POSITIVE,'LineStyle','-','LineWidth',LINE_WIDTH*5);
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
        if isempty(labels), continue; end
        maxx = max(scan.Vertices(1,:));
        minx = min(scan.Vertices(1,:));
        pointx = maxx(1)-(3/4)*(maxx(1)-minx(1));
        pointz = max(scan.Vertices(3,:));pointz = pointz(1);
        maxy = max(scan.Vertices(2,:));
        miny = min(scan.Vertices(2,:));
        pointy = (maxy(1)+miny(1))/2;
        %pointx = min(scan.Vertices(1,:));pointx = pointx(1);
        %[~,maxpoint] = max(scan.Vertices(3,:));
        %maxpoint = scan.Vertices(:,maxpoint);
        %maxpoint(3) = maxpoint(3)+0.5;
        %text(maxpoint(1),maxpoint(2),maxpoint(3),labels{n},'Color','k','FontSize',20,'Parent',plotax);
        text(pointx,pointy,pointz,labels{n},'Color','k','FontSize',10,'Parent',plotax);
    end
    view(plotax,0,90);
    %plotax.Visible = 'off';
    plotax.Color = [1 1 1];
    plotax.XColor = [1 1 1];plotax.YColor = [1 1 1];plotax.ZColor = [1 1 1];
    colormap(plotax,cm);colorbar('peer',plotax);
    set(plotax,'clim',range);
    % Set the axis appropriately
    max_space = MARKER_MAX_SIZE/2;
    xlim([min(point_list(1,:))-max_space,max(point_list(1,:))+max_space]);
    ylim([min(point_list(2,:))-max_space,max(point_list(2,:))+max_space]);
    switch v.Type
           case 'axes'
               axis(plotax,'square');
           case 'viewer3DObj'
               axis(plotax,'equal');
           case 'figure'
               axis(plotax,'square');
    end
    %hold(plotax,'off');
    if ~isempty(titlestr), title(plotax,titlestr,'Color',[0 0 0]);end
    if ~(light),return; end
    switch v.Type
           case 'axes'
               camlight headlight;
           case 'viewer3DObj'
               %camlight headlight;
               v.SceneLightVisible = true;
               v.SceneLightLinked = true;
           case 'figure'
               camlight headlight;
    end
end

