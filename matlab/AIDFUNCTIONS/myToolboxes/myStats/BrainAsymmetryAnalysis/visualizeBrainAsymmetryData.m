function f = visualizeBrainAsymmetryData(brainSurface, data, nSamplesPerPick, titlenames)
    if isstring(data)
        load(data,'data' ,'nSamplesPerPick', 'titlenames', 'brainSurface');
    end
    VertexValues = data;
    

    arrange = [3 8];
    counter = 0;
    clim = [];
    rend = brainSurface;
    cmap = parula(256);

    if size(VertexValues{1,1},1) > 1
        for i= 1:3
            f = showSignificanceEvolution(brainSurface, nSamplesPerPick, VertexValues{i,4},  strcat(titlenames{i,1}," ",titlenames{i,4}));
        end
        return
    end
    dmap = parula(length(thresholds));

    f = figure;f.Position = [95  98  2192  1106];f.Color = [1 1 1];%
    thresholdsTicks = thresholds;
    for i=1:nValues
        for j=1:4
            if j==4, map=dmap; else, map=cmap; end
            %i=1;
            counter = counter+1;
            fout.ax1{i,j} = subplot(arrange(1),arrange(2),counter,'Parent',f);
            colormap(fout.ax1{i,j},map);

            if ~isempty(clim), set(fout.ax1{i,j},'clim',clim);end
            renderBrainSurface(rend,VertexValues{i,j},fout.ax1{i,j});
            colorbar(fout.ax1{i,j},'SouthOutside');
            if j<4
                set(fout.ax1{i,j},'clim',[0 max(VertexValues{i,j})]);
            end
            if j==4,set(fout.ax1{i,j},'clim',[0 length(thresholds)]); cb=findall(gcf,'type','ColorBar');cb(1).Ticks=1:length(thresholdsTicks);
                cb(1).TickLabels=thresholdsTicks; cb(1).Label.String = 'p-value'; end
            view(fout.ax1{i,j},rend.viewval(1),0);
            light = camlight(fout.ax1{i,j},'headlight');
            set(light,'Position',get(fout.ax1{i,j},'CameraPosition'));
            drawnow;
            title(fout.ax1{i,j},titlenames{i,j})
            counter = counter+1;
            fout.ax2{i,j} = subplot(arrange(1),arrange(2),counter,'Parent',f);
            renderBrainSurface(rend,VertexValues{i,j},fout.ax2{i,j});
            view(fout.ax2{i,j},-1*rend.viewval(1),0);
            colorbar(fout.ax2{i,j},'SouthOutside');
            if j<4
                set(fout.ax2{i,j},'clim',[0 max(VertexValues{i,j})]);
            end
            colormap(fout.ax2{i,j},map);
            if j==4,set(fout.ax2{i,j},'clim',[0 length(thresholds)]); cb=findall(gcf,'type','ColorBar');cb(1).Ticks=1:length(thresholdsTicks);
                cb(1).TickLabels=thresholdsTicks; cb(1).Label.String = 'p-value'; end
            light = camlight(fout.ax2{i,j},'headlight');
            set(light,'Position',get(fout.ax2{i,j},'CameraPosition'));
            drawnow;
            if ~isempty(clim), set(fout.ax2{i,j},'clim',clim);end
        end
    end

    end

