function [out,compactnessPerLM,AvgLM,DATA,time] = compactnessMapping(DATA,RMapper,NRMapper,Template,observer,display)
    tic;
    % run the mapping on all the DATA
    N = length(DATA);
    [path,ID] = setupParForProgress(N);
    parfor i=1:N
        if ~isempty(RMapper), forRMapper = clone(RMapper); else, forRMapper = []; end
        if ~isempty(NRMapper), forNRMapper = clone(NRMapper); else, forNRMapper = []; end
        out = mapFace(DATA{i}.Shape,Template,forRMapper,forNRMapper,false);
        DATA{i}.MappedShape = clone(out);
        parfor_progress;
    end    
    closeParForProgress(path,ID);
    % extract landmarks
    for i=1:1:N
        [bar,index] = cart2baryKNN(DATA{i}.MappedShape.Vertices,DATA{i}.AvgLM{observer}.Vertices);
        MappedLM = clone(DATA{i}.AvgLM{observer});
        [MappedLM.Vertices] = bary2cartKNN(Template.Vertices,index,bar);        
        DATA{i}.BAR = bar;
        DATA{i}.index = index;
        DATA{i}.MappedLM = MappedLM;
    end
    if display, v = viewer(clone(Template)); end
    % computing average landmarks
    AvgLM = clone(DATA{1}.AvgLM);
    AvgLM.SingleColor = [0 1 0];
    LM = zeros(AvgLM.nVertices,3,N);
    for i=1:1:N
        LM(:,:,i) = DATA{i}.MappedLM.Vertices;
    end
    AvgLM.Vertices = mean(LM,3);
    if display
       AvgLM.VertexSize =40;
       viewer(AvgLM,v);
       for i=1:1:N
           DATA{i}.MappedLM.VertexSize = 10;
           viewer(DATA{i}.MappedLM,v);
       end
    end
    % projecting onto the surface
    [bar,index] = cart2baryKNN(Template.Vertices,AvgLM.Vertices);
    [AvgLM.Vertices] = bary2cartKNN(Template.Vertices,index,bar);
    % compute compactness per landmark
    nLM = AvgLM.nVertices;
    compactnessPerLM = zeros(1,nLM);
    for l=1:1:nLM
        tmp = squeeze(LM(l,:,:))';
        dist = pdist2(AvgLM.Vertices(l,:),tmp);
        compactnessPerLM(l) = mean(dist);
    end
    if display
       AvgLM.VertexValue = compactnessPerLM;
       AvgLM.ColorMode = 'Indexed';
       set(gca,'Clim',[0 3]);
       v.ColorbarVisible = 'on';
    end   
    out = sum(compactnessPerLM);
    time = toc;
end