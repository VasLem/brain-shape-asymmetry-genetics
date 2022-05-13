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
arrange = [size(VertexValues,1) size(VertexValues,2)];

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
            fout.ax{i,j} = subplot(arrange(1),arrange(2),counter,'Parent',f);
            view(fout.ax{i,j},-90,0);
            shape = clone(brainSurface);
            shape.ColorMode = 'indexed';
            shape.VertexValue = sVertexValues;
            shape.Material = 'dull';

            if ~isempty(clim), set(fout.ax{i,j},'clim',clim);end
            showPaintedDoubleFace(fout, shape, nan, nan, fout.ax{i,j})
            colormap(fout.ax{i,j}, map)
            colorbar(fout.ax{i,j},'SouthOutside');
            if j<4
                set(fout.ax{i,j},'clim',[0 max(sVertexValues)]);
            end
            if j==4
                set(fout.ax{i,j},'clim',[0 length(thresholds)]); 
                cb=findall(f,'type','ColorBar');
                cb(1).Ticks=1:length(thresholdsTicks);
                cb(1).TickLabels=thresholdsTicks;
                cb(1).Label.String = 'p-value';
                cb(1).Label.Interpreter = 'latex';
            end
            drawnow;
            fout.ax{i,j}.Visible = 'on';
            axis(fout.ax{i,j},'image');
            axis(fout.ax{i,j},'off');
            title(fout.ax{i,j},titlenames{i,j},'FontUnit','normalized', 'FontSize',0.2);
            
            
            figName = saveName;
            if nPicks > 1, figName = [saveName '_' num2str(k)]; end
               
            saveas(f, [figName '.fig']);
            saveas(f, [figName '.svg']);
        end
    end
    
end

