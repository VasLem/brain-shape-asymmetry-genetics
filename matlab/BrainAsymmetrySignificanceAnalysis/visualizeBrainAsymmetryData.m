function f = visualizeBrainAsymmetryData(data,saveName)
if isstring(data)
    load(data,'data');
end
VertexValues = data.values;
brainSurface = data.brainSurface;
titlenames = data.titleNames;
try
    thresholds = data.thresholds;
    nSamplesPerPick = data.nSamplesPerPick;
catch
end
arrange = [size(VertexValues,1) 2 * size(VertexValues,2)];

clim = [];
rend = brainSurface;
cmap = parula(256);
dmap = 0;
if size(VertexValues{1,1},1) > 1
    for i= 1:3
        showSignificanceEvolution(brainSurface, nSamplesPerPick, VertexValues{i,4},  strcat(titlenames{i,1}," ",titlenames{i,4}));
    end
end
if size(VertexValues,2) > 1
    dmap = parula(length(thresholds));
    thresholdsTicks = thresholds;
end

nPicks = size(VertexValues{1,1},1);
for k=1:nPicks
    counter = 0;
    try
        f = figure('Name', ['Picked:' num2str(nSamplesPerPick(k))]);
    catch
        f = figure;
    end
    f.Position = [95  98  2192  1106];f.Color = [1 1 1];%
    for i=1:size(VertexValues,1)
        for j=1:size(VertexValues,2)
            if j==4, map=dmap; else, map=cmap; end
            %i=1;
            sVertexValues = VertexValues{i,j}(k,:);
            counter = counter+1;
            fout.ax1{i,j} = subplot(arrange(1),arrange(2),counter,'Parent',f);
            colormap(fout.ax1{i,j},map);
            
            if ~isempty(clim), set(fout.ax1{i,j},'clim',clim);end
            renderBrainSurface(rend,sVertexValues,fout.ax1{i,j});
            maxlims = max(rend.Vertices);
            minlims = min(rend.Vertices);
            colorbar(fout.ax1{i,j},'SouthOutside');
            if j<4
                set(fout.ax1{i,j},'clim',[0 max(sVertexValues)]);
            end
            if j==4,set(fout.ax1{i,j},'clim',[0 length(thresholds)]); cb=findall(gcf,'type','ColorBar');cb(1).Ticks=1:length(thresholdsTicks);
                cb(1).TickLabels=thresholdsTicks; cb(1).Label.String = 'p-value'; end
            view(fout.ax1{i,j},-90,0);
            xlim(fout.ax1{i,j}, [0, 1.3 * (maxlims(1) + minlims(1)) / 3])
            light = camlight(fout.ax1{i,j},'headlight');
            set(light,'Position',get(fout.ax1{i,j},'CameraPosition'));
            drawnow;
            title(fout.ax1{i,j},titlenames{i,j})
            counter = counter+1;
            fout.ax2{i,j} = subplot(arrange(1),arrange(2),counter,'Parent',f);
            renderBrainSurface(rend,sVertexValues,fout.ax2{i,j});

            view(fout.ax2{i,j},90,0);
            colorbar(fout.ax2{i,j},'SouthOutside');
            if j<4
                set(fout.ax2{i,j},'clim',[0 max(sVertexValues)]);
            end
            colormap(fout.ax2{i,j},map);
            if j==4,set(fout.ax2{i,j},'clim',[0 length(thresholds)]); cb=findall(gcf,'type','ColorBar');cb(1).Ticks=1:length(thresholdsTicks);
                cb(1).TickLabels=thresholdsTicks; cb(1).Label.String = 'p-value'; end
            light = camlight(fout.ax2{i,j},'headlight');
            set(light,'Position',get(fout.ax2{i,j},'CameraPosition'));
            drawnow;
            if ~isempty(clim), set(fout.ax2{i,j},'clim',clim);end
            figName = saveName;
            if nPicks > 1, figName = [saveName '_' num2str(k)]; end
               
            saveas(f, [figName '.fig']);
            saveas(f, [figName '.svg']);
        end
    end
    
end

