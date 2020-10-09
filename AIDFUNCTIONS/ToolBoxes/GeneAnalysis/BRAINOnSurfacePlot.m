function [v] = BRAINOnSurfacePlot(Render,segmentvalues,level)
    VAL = nan*zeros(size(Render.ULABELS));
    counter = 0;
    for i=1:Render.UHI.nLC
       if Render.UMASK(i)==0, continue;end
       counter = counter+1;
       % i=3
       [l,c] = Ind2LC(Render.UHI,i);
       index = find(Render.ULABELS(l,:)==c);
       VAL(l,index) = segmentvalues(counter);
    end
    if level>size(VAL,1), level = size(VAL,1); end
    VAL = VAL(1:level,:);
    
    values = VAL(end,:);
    counter = 1;
    while ~isempty(find(isnan(values)))&&counter<level
        index = find(isnan(values));
        values(index) = VAL(end-counter,index);
        counter = counter+1;
    end
    
    RefScan = clone(Render.RefScan);
    RefScan.VertexValue = values;
    RefScan.ColorMode = "Indexed";
    RefScan.Material = 'Dull';
    v = viewer(RefScan);
    v.Figure.Color = [1 1 1];
    v.AxesXColor = [0 0 0];
    v.AxesYColor = [0 0 0];
    v.AxesZColor = [0 0 0];
    v.Colorbar.Color = [0 0 0];
    view(Render.viewval(1),0);
    v.SceneLightVisible = true;
    v.SceneLightLinked = true;
    RefScan.PatchHandle.FaceColor = 'flat';
end