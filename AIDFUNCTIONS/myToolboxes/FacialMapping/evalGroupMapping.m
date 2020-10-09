function [out,compactnessPerLM,AvgLM,LM] = evalGroupMapping(DATA,Template,observer,display,LandmarkNames)
    N = length(DATA);
    % extract landmarks
    for i=1:1:N
        %i=1;
        [bar,index] = cart2baryKNN(DATA{i}.MappedShape.Vertices,DATA{i}.AvgLM{observer}.Vertices);
        MappedLM = clone(DATA{i}.AvgLM{observer});
        [MappedLM.Vertices] = bary2cartKNN(Template.Vertices,index,bar);
        DATA{i}.LMOnTemplate = MappedLM;
    end
    if display, v = viewer(clone(Template)); end
    % computing average landmarks
    AvgLM = clone(DATA{1}.AvgLM{observer});
    AvgLM.SingleColor = [0 1 0];
    LM = zeros(AvgLM.nVertices,3,N);
    for i=1:1:N
        LM(:,:,i) = DATA{i}.LMOnTemplate.Vertices;
    end
    AvgLM.Vertices = mean(LM,3);
    if display
       AvgLM.VertexSize =30;
       viewer(AvgLM,v);
       for i=1:1:N
           DATA{i}.LMOnTemplate.VertexSize = 10;
           DATA{i}.LMOnTemplate.SingleColor = [1 0 0];
           viewer(DATA{i}.LMOnTemplate,v);
       end
    end
    % projecting onto the surface
    [bar,index] = cart2baryKNN(Template.Vertices,AvgLM.Vertices);
    [AvgLM.Vertices] = bary2cartKNN(Template.Vertices,index,bar);
    % compute compactness per landmark
    nLM = AvgLM.nVertices;
    compactnessPerLM = zeros(1,nLM);
    dist = zeros(N,nLM);
    for l=1:1:nLM
        tmp = squeeze(LM(l,:,:))';
        dist(:,l) = pdist2(AvgLM.Vertices(l,:),tmp);
        compactnessPerLM(l) = mean(dist(:,l));
    end
    if display
       AvgLM.VertexValue = compactnessPerLM;
       AvgLM.ColorMode = 'Indexed';
       set(gca,'Clim',[0 3]);
       v.ColorbarVisible = 'on';
       f = figure;boxplot(dist);grid on;set(gca,'ylim',[0 10]);
       f.CurrentAxes.XTickLabel = LandmarkNames(:,1);
       f.CurrentAxes.XTickLabelRotation = 90;
       f.Color = [1 1 1];
       ylabel('(mm)');
       title('Landmark Compactness after templating');
    end   
    out = sum(compactnessPerLM);
end