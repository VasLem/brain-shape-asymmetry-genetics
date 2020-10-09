function peer = plotHierarchicalLabelsDGors(peer,rend,SegmentVal, caxis_min_max, clmap, crit)
% This function plots a circular dendrogram from the results calculated within a hierarchical segmentation.
%
%   use:
%           D_plotHierarchicalLabels_CircularDendrogram(labels);
%           D_plotHierarchicalLabels_CircularDendrogram(labels, caxis_min_max, clmap, minWhite, circleSquare, TitleColorbar);
%   inputs:
%           labels =            nLevels x 1 cell: {[x];[x x];[x x x x];...}
%           caxis_min_max =     [minValue of colorbar, maxValue of colorbar]: set a special range for the colorbar
%           clmap :             default 'jet'
%           minWhite:           [] (default), true (= below minimum is collored grey) (works only with the correct caxis_min_max inputs)
%           circleSquare:       [] (default): nodes are plotted as circles
%                               nLevels x 1 cell: {[x];[x x];[x x x x];...}: x = 0 (circles), x = 1 (squares)
%% TURN VALUES INGO LABELS
    MASK = rend.UMASK;
    HI = rend.UHI;
    % I need to transform my values into labels
    nLevels = HI.nL;
    values = nan*zeros(1,HI.nLC);
    values(find(MASK)) = SegmentVal;
    labels = cell(nLevels,1);
    for i=1:HI.nLC
       [l,c] = Ind2LC(HI,i);
       labels{l}(c) = values(i);
    end
    if isempty(crit)
        circleSquare = [];
    else
        circleSquare = cell(nLevels,1);
        for i=1:HI.nLC
           [l,c] = Ind2LC(HI,i);
           circleSquare{l}(c) = values(i)>=crit;
        end
    end
    %%
    if exist('clmap','var')==1 && isempty(clmap)==0
        %Nodes_color = eval([clmap,'(10000);']);
        Nodes_color = clmap;
    else
        %Color_Nodes = 'jet';
        Nodes_color = colormap(peer,'parula');
    end
    nColors = size(clmap,1);
%     if exist('minWhite','var')==1 && isempty(minWhite)==0 && minWhite == true
%         Nodes_color(1,:) = 0.6;
%     end
    if exist('circleSquare','var')==0 || isempty(circleSquare)
        circleSquare = cell(nLevels,1);
        for level = 1:nLevels
            circleSquare{level} = zeros(1,2^(level-1));
        end
    end
    nClusters = sum(2.^[0:nLevels-1]);       
    id2LevelCluster = zeros(nClusters,2);
    count = 0;
    for level = 1:nLevels
        nCl = 2^(level-1);
        id2LevelCluster(count+1:count+nCl,:) = [level*ones(nCl,1), [1:nCl]'];
        count = count + nCl;
    end
    %% Set the variables
    color_HierCircles = [0.6 0.6 0.6];
    color_HierLines = [0.6 0.6 0.6];
    lineWidth_HierCircles = 0.5;
    lineWidth_HierLines = 0.5;
    %% Plot the outher circle
    %if isempty(peer); figure; peer = gca;end
    %axes(peer);
    %hold on
    hold(peer,'on');
    colormap(peer,'jet')
    RADIUS = 1;
    rectangle(peer,'Position',[-RADIUS, -RADIUS, 2*RADIUS, 2*RADIUS],...
        'Curvature', [1 1],'EdgeColor',color_HierCircles,'LineWidth',lineWidth_HierCircles,'LineStyle',':','FaceColor','none');
    %% Calculate the dendrogram range
    if exist('caxis_min_max','var')==0 || isempty(caxis_min_max)==1,
        caxis_min_max = [min(cell2mat(labels')) max(cell2mat(labels'))];
    end
    %% Plot the hierarchical structure
    nodes_levels = cell(nLevels,1);
    nodes_levels{1} = [0 0];
    for level = 2:nLevels,
        % draw the child-circle
        RADIUS_level = (level-1)*RADIUS/(nLevels-1);
        rectangle(peer,'Position',[-RADIUS_level, -RADIUS_level, 2*RADIUS_level, 2*RADIUS_level],...
            'Curvature', [1 1],'EdgeColor',color_HierCircles,'LineWidth',lineWidth_HierCircles,'LineStyle',':','FaceColor','none');
        nCl = 2^(level-1);
        degree_space = 360/nCl;
        %degree_serie = -90+[degree_space/2:degree_space:360-0.5*degree_space]';
        degree_serie = -90+[360-0.5*degree_space:-degree_space:degree_space/2]';
        nodes_levels{level} = RADIUS_level*[cosd(degree_serie) -sind(degree_serie)];
        % draw the connection to the parent
        for cl = 1:nCl,
            [ind] = LC2Ind(HI,level,cl);
            if MASK(ind)==1
                cl_parent = ceil(cl/2);
                plot(peer,[nodes_levels{level-1}(cl_parent,1) nodes_levels{level}(cl,1)],[nodes_levels{level-1}(cl_parent,2) nodes_levels{level}(cl,2)],'-','Color',color_HierLines,'LineWidth',lineWidth_HierLines)
            end
        end
    end
    %% plot the collored Nodes
    for level = 1:nLevels;
        for cl = 1:2^(level-1),
            score = labels{level}(cl);
            if ~isnan(score);
                if score > caxis_min_max(1),
                    if score < caxis_min_max(2)
                        score = (score-caxis_min_max(1))/(caxis_min_max(2)-caxis_min_max(1));
                        if round(nColors*score)==0; score = 1/nColors; end
                        %nodekleur = Nodes_color(round(nColors*score),:);
                        %nodekleur = Nodes_color(round(nColors*score),:);
                        nodekleur = Nodes_color(ceil(nColors*score),:);
                    else
                        nodekleur = Nodes_color(end,:);
                    end
                else
                    nodekleur = Nodes_color(1,:);
                end
                
                max_NodeRadius = 0.062;
                if level > 1;
                    dd = 0.5*pdist2(nodes_levels{level}(1,:),nodes_levels{level}(2,:));
                    if dd < max_NodeRadius,
                        NodeRadius = 0.85*dd;
                    else
                        NodeRadius = max_NodeRadius;
                    end
                else
                     NodeRadius = max_NodeRadius;
                end
                    
                if ~isequal(nodekleur,[0.6 0.6 0.6])
                    if circleSquare{level}(cl) == 0;
                        rectangle(peer,'Position',[nodes_levels{level}(cl,1)-NodeRadius, nodes_levels{level}(cl,2)-NodeRadius, 2*NodeRadius, 2*NodeRadius],...
                            'Curvature', [1 1],'EdgeColor','none','FaceColor',nodekleur);
                    else
                        NodeRadius = 0.75*NodeRadius;
                        rectangle(peer,'Position',[nodes_levels{level}(cl,1)-NodeRadius, nodes_levels{level}(cl,2)-NodeRadius, 2*NodeRadius, 2*NodeRadius],...
                            'Curvature', [0 0],'EdgeColor','none','FaceColor',nodekleur);
                    end
                else
                    NodeRadius = 0.75*NodeRadius;
                        rectangle(peer,'Position',[nodes_levels{level}(cl,1)-NodeRadius, nodes_levels{level}(cl,2)-NodeRadius, 2*NodeRadius, 2*NodeRadius],...
                            'Curvature', [0.6 0.6],'EdgeColor','none','FaceColor',nodekleur);
                end
            end
        end
    end
    %axis(peer,'square');
    axis(peer,'image');
    %axis(peer,[-1.15 1.15 -1.15 1.15]);
    axis(peer,'off');
    caxis(peer,caxis_min_max);
    colormap(peer,clmap);
end

