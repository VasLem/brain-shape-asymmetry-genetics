function f = D_plotConnectome( Connectome, Box_values, Range_Connectome, Range_Box, strarray, Title2Plot)
% Plot a connectome.
%           D_plotConnectome(Connectome);
% Inputs: 
%           Connectome: nClusters x nClusters-matrix with the weights of the connectome-connections
%           Box_values: nClusters-array (optional) with the intracluster weights
%           Range_Connectome: plot range for the connectome-connections (optional)
%           Range_Box: plot range for the intracluster weights (optional)
%           Hierarchical_Level: level of the hierarchical tree. If given the hierarchical tree will be plotted under the connectome structure
%           Title2Plot: Title for the figure (optional)
%
% To plot a disconnectome ranging from a to b (a>b):
%           D_plotConnectome(1-Connectome,[],[1-a 1-b])
%           cb = colorbar('location','eastoutside','Ticks',[0:0.1:1],'TickLabels',{[linspace(a,b,11)]});
%           title(cb,'connections')

    nClusters = size(Connectome,1);
    if exist('Hierarchical_Level','var')==1 && isempty(Hierarchical_Level)==0,
        doHier = true;
        if nClusters ~= 2^(Hierarchical_Level-1),
            fprintf('D_plotConnectome: The number of clusters do not match with the level.\n');
        end
    else
        doHier = false;
    end

    %% Set the variables
    max_LineWidth = 3;
    max_BoxRadius = 0.1;
    color_BoxDefault = [0 0.6 1];
    color_HierCircles = [0.5 0.5 0.5];
    color_HierLines = [0.5 0.5 0.5];
    lineWidth_HierCircles = 1;
    lineWidth_HierLines = 2;
    linecolors_jet = parula(1000);
    linecolors_summer = summer(1000);
    %linecolors_summer = lines(1000);
    linecolors_lines = lines;
    index = find(ismember(linecolors_lines,[0.4660    0.6740    0.1880],'rows'));
    for i=1:length(index)
        linecolors_lines(index(i),:) = [0.4 0.4 0.4];
    end

    %% Plot the outher circle
    f = figure();
    hold on
    colormap ('jet')
    RADIUS = 1;
    rectangle('Position',[-RADIUS, -RADIUS, 2*RADIUS, 2*RADIUS],...
    'Curvature', [1 1],'EdgeColor',color_HierCircles,'LineWidth',lineWidth_HierCircles,'LineStyle',':','FaceColor','none');

    %% Plot the hierarchical structure
    if doHier==1,
        nodes_levels = cell(Hierarchical_Level,1);
        nodes_levels{1} = [0 0];
        for level = 2:Hierarchical_Level,
            % draw the child-circle
            RADIUS_level = (level-1)*RADIUS/(Hierarchical_Level-1);
            rectangle('Position',[-RADIUS_level, -RADIUS_level, 2*RADIUS_level, 2*RADIUS_level],...
                'Curvature', [1 1],'EdgeColor',color_HierCircles,'LineWidth',lineWidth_HierCircles,'LineStyle',':','FaceColor','none');
            
            nCl = 2^(level-1);
            degree_space = 360/nCl;
            degree_serie = -90+[degree_space/2:degree_space:360-0.5*degree_space]';
            nodes_levels{level} = RADIUS_level*[cosd(degree_serie) -sind(degree_serie)];
            % draw the connection to the parent
            for cl = 1:nCl,
                cl_parent = ceil(cl/2);
                plot([nodes_levels{level-1}(cl_parent,1) nodes_levels{level}(cl,1)],[nodes_levels{level-1}(cl_parent,2) nodes_levels{level}(cl,2)],'-','Color',color_HierLines,'LineWidth',lineWidth_HierLines)
            end
        end
        nodes_pos = nodes_levels{level};    clear nodes_levels;
        
    else
        degree_space = 360/nClusters;
        degree_serie = -90+[degree_space/2:degree_space:360-0.5*degree_space]';
        %nodes_pos = RADIUS*[cosd(nodes_degree) -sind(nodes_degree)];
        nodes_pos = RADIUS*[cosd(degree_serie) -sind(degree_serie)];
        nodes_pos2 = RADIUS*1.025*[cosd(degree_serie) -sind(degree_serie)];
    end
    
    %% Plot the connectome lines
    if exist('Range_Connectome','var')==0 || isempty(Range_Connectome)==1,
        if (exist('Box_values','var')==1 && isempty(Box_values)==0) && (exist('Range_Box','var')==0 || isempty(Range_Box)==1),
            Range_Connectome = [min([Connectome(:);Box_values(:)]) max([Connectome(:);Box_values(:)])];
        else
            Range_Connectome = [min(Connectome(:)) max(Connectome(:))];
        end
    end
    for cl1 = 1:nClusters,
        %cl1 = 1;
        for cl2 = cl1+1:nClusters,
            %cl2 = 2;
            score = Connectome(cl1,cl2);
            if score >= Range_Connectome(1),
                if score <= Range_Connectome(2),
                    score = (score-Range_Connectome(1))/(Range_Connectome(2)-Range_Connectome(1));
                    %lijndikte = max(0.5, score*max_LineWidth);
                    lijndikte = max(0.5, max_LineWidth);
                    
                    if round(1000*score)==0; score = 1/1000; end
                    lijnkleur = [linecolors_jet(round(1000*score),:) score];
                    
                    %plot([nodes_pos(cl1,1) nodes_pos(cl2,1)],[nodes_pos(cl1,2) nodes_pos(cl2,2)],'-','linewidth',lijndikte,'color',lijnkleur)
                    % determine the centre and radius of connection line
                    
                    node1 = nodes_pos(cl1,:);
                    node2 = nodes_pos(cl2,:);
                    % determine the centre and radius of connection line
                    if cl2-cl1 ~= nClusters/2,
                        x = [-node1(2) node2(2);node1(1) -node2(1)]\[-node1(1)+node2(1);-node1(2)+node2(2)];
                        node3 = node1'+x(1)*[-node1(2);node1(1)];
                        node3 = node3';
                        radius = pdist2(node1,node3);        
                        
                        % plot enkel tussen de twee nodes
                        hoek1  = atan2d(node1(2)-node3(2),node1(1)-node3(1));   
                        hoek2  = atan2d(node2(2)-node3(2),node2(1)-node3(1)); 
                        if abs(hoek1-hoek2)>180,
                            if hoek1<0, hoek1 = hoek1+360; end
                            if hoek2<0, hoek2 = hoek2+360; end
                        end
%                         if sign(hoek1)~=sign(hoek2);
%                             if hoek1<0, hoek1 = hoek1+360; end
%                             if hoek2<0, hoek2 = hoek2+360; end
%                         end
                        hoek_connection = linspace(hoek1,hoek2,100);
                        x_connection = node3(1) + radius*cosd(hoek_connection);
                        y_connection = node3(2) + radius*sind(hoek_connection);
                        plot(x_connection,y_connection,'-','linewidth',lijndikte,'color',lijnkleur)
                        
                    else
                        plot([nodes_pos(cl1,1) nodes_pos(cl2,1)],[nodes_pos(cl1,2) nodes_pos(cl2,2)],'-','linewidth',lijndikte,'color',lijnkleur)
                    end
                else
                    % do not plot to strong connections
                end
            else
                % do not plot to small connections
            end
        end
    end
    
    %% determine of the box radius needs to be adjusted
    dd = 0.5*pdist2(nodes_pos(1,:),nodes_pos(2,:));
    if dd < max_BoxRadius,
    max_BoxRadius = 0.85*dd;
    end
    
    %% plot the collored boxes
    if exist('Box_values','var')==1 && isempty(Box_values)==0,
        if exist('Range_Box','var')==0 || isempty(Range_Box)==1,
            Range_Box = Range_Connectome;
        end
        for cl = 1:nClusters,
            score = Box_values(cl);
            if score >= Range_Box(1),
                if score <= Range_Box(2)
                    score = (score-Range_Box(1))/(Range_Box(2)-Range_Box(1));
                    if round(1000*score)==0; score = 1/1000; end
                    boxkleur = linecolors_summer(round(1000*score),:);
                else
                    boxkleur = linecolors_summer(end,:);
                end
            else
                boxkleur = linecolors_summer(1,:);
            end
            
            boxkleur = linecolors_lines(Box_values(cl),:);
            
            rectangle('Position',[nodes_pos(cl,1)-max_BoxRadius, nodes_pos(cl,2)-max_BoxRadius, 2*max_BoxRadius, 2*max_BoxRadius],...
                'Curvature', [1 1],'EdgeColor','none','FaceColor',boxkleur);
        end
    else
        for cl = 1:nClusters,
            rectangle('Position',[nodes_pos(cl,1)-max_BoxRadius, nodes_pos(cl,2)-max_BoxRadius, 2*max_BoxRadius, 2*max_BoxRadius],...
                'Curvature', [1 1],'EdgeColor','none','FaceColor',color_BoxDefault);
        end
    end
    
    axis square
    axis ([-1.15 1.15 -1.15 1.15])
    axis off
    
    if exist('strarray','var')==1
       for i=1:length(strarray)
           %i=3;
           if ~isempty(strarray{i})
               text(nodes_pos2(i,1),nodes_pos2(i,2),strarray{i},'HorizontalAlignment','left','FontSize',6,'Rotation',-1*degree_serie(i));
           end
       end
    end
    
    ticks = 0:0.25:1;
    subval = round((Range_Connectome(2)-Range_Connectome(1))/4);
    %ticklabels = [Range_Connectome(1) round((Range_Connectome(2)-Range_Connectome(1))/2) Range_Connectome(2)];
    %ticklabels = [Range_Connectome(1) Range_Connectome(1)+subval Range_Connectome(1)+2*subval Range_Connectome(1)+3*subval Range_Connectome(2)];
    ticklabels = round([1 Range_Connectome(1)+subval Range_Connectome(1)+2*subval Range_Connectome(1)+3*subval Range_Connectome(2)]);
    %cb = colorbar('location','eastoutside','Limits',[Range_Connectome(1) Range_Connectome(2)], 'Ticks',[Range_Connectome(1):1:Range_Connectome(2)],'TickLabels',[Range_Connectome(1):1:Range_Connectome(2)]);
    cb = colorbar('location','eastoutside', 'Ticks',ticks,'TickLabels',ticklabels);
    colormap(cb,'parula');
    title(cb,'Modules Connected')
    if exist('Box_values','var')==1 && isempty(Box_values)==0
        ticks = 0:0.0435:1;
        ticklabels = cell(1,23);
        for i=1:22
            ticklabels{i} = num2str(i);
        end
        ticklabels{end} = 'X';
        cb = colorbar('southoutside','Ticks',ticks,'TickLabels',ticklabels);
        colormap(cb,linecolors_lines(1:23,:));
        title(cb,'Chromosome')
    end
    
    set(gcf,'position',[130   131   967   888])
    if exist('Title2Plot','var')==1 && isempty(Title2Plot)==0
    title(Title2Plot);
    end
    set(gca,'fontsize',14)



end

