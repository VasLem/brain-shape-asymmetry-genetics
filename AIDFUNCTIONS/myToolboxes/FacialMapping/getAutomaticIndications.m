function [out,compactnessPerLM,AvgLM,LM,DATA] = getAutomaticIndications(DATA,Template,display,LandmarkNames)
    N = length(DATA);
    % FIRST: TRANSFORM ALL THE LANDMARKS TO TEMPLATE SPACE
    for i=1:1:N
        %i=1;
        DATA{i}.AvgLMOnTemplate = cell(1,3);
        for o=1:3% for all three observers
           [bar,index] = cart2baryKNN(DATA{i}.MappedShape.Vertices,DATA{i}.AvgLM{o}.Vertices);
           MappedLM = clone(DATA{i}.AvgLM{o});
           [MappedLM.Vertices] = bary2cartKNN(Template.Vertices,index,bar);        
           DATA{i}.AvgLMOnTemplate{o} = MappedLM;
        end
    end
    if display, v = viewer(clone(Template)); end
    % computing average landmarks
    LM = zeros(DATA{1}.AvgLM{1}.nVertices,3,N,3);
    for i=1:1:N
        for o=1:1:3
            LM(:,:,i,o) = DATA{i}.AvgLMOnTemplate{o}.Vertices;
        end
    end
    AvgLM = cell(1,3);
    AvgBar = cell(1,3);
    AvgIndex = cell(1,3);
    for o=1:1:3
       AvgLM{o} = clone(DATA{1}.AvgLM{1});
       AvgLM{o}.Vertices = mean(squeeze(LM(:,:,:,o)),3);
       c = zeros(1,3);c(o) = 1;
       AvgLM{o}.SingleColor = c; 
       AvgLM{o}.VertexSize =30;
       
       % projecting onto the surface
       [AvgBar{o},AvgIndex{o}] = cart2baryKNN(Template.Vertices,AvgLM{o}.Vertices);
       [AvgLM{o}.Vertices] = bary2cartKNN(Template.Vertices,AvgIndex{o},AvgBar{o});
       if display, viewer(AvgLM{o},v);end 
    end  
    % compute compactness per landmark
    nLM = AvgLM{1}.nVertices;
    compactnessPerLM = zeros(nLM,3);
    for l=1:1:nLM
        for o=1:1:3
            tmp = squeeze(LM(l,:,:,o))';
            dist = pdist2(AvgLM{o}.Vertices(l,:),tmp);
            compactnessPerLM(l,o) = mean(dist);
        end
    end
    out = sum(compactnessPerLM,1);
    % SECOND: TRANFORM AVERAGES TO THE MESHES
    for i=1:1:N
        %i=1;
        DATA{i}.AvgLMOnFace = cell(1,3);
        for o=1:3% for all three observers
           MappedLM = clone(DATA{i}.AvgLM{o});
           [MappedLM.Vertices] = bary2cartKNN(DATA{i}.MappedShape.Vertices,AvgIndex{o},AvgBar{o});
           % maybe we need to project further
           [bar,index] = cart2baryKNN(DATA{i}.Shape.Vertices,MappedLM.Vertices);
           [MappedLM.Vertices] = bary2cartKNN(DATA{i}.Shape.Vertices,index,bar);
           DATA{i}.AvgLMOnFace{o} = MappedLM;
        end
    end
%     for i=1:1:N
%        v = viewer(DATA{i}.Shape);
%        for o=1:1:3
%            viewer(DATA{i}.AvgLMOnFace{o},v);
%        end
%        waitfor(v);
%     end
    % THIRD:MAKE MEASERMENTS
    dist = zeros(AvgLM{1}.nVertices,N,3,3);
    %i=1;o=1;oo=1;
    for i=1:N
       for o=1:1:3
           for oo=1:1:3
               dist(:,i,o,oo) = sqrt(sum((DATA{i}.AvgLM{o}.Vertices-DATA{i}.AvgLMOnFace{oo}.Vertices).^2,2));
           end
       end
    end
    
    
%     nLM = AvgLM{1}.nVertices;
%     DIST = zeros(N*3*3,nLM);
%     for l=1:1:nLM
%         tmp = squeeze(dist(l,:,:,:));
%         DIST(:,l) = tmp(:);
%     end
%     index = find(DIST(:,1));
%     DIST = DIST(index,:);
%     
%     f = figure;boxplot(DIST);grid on;set(gca,'ylim',[0 10]);
    
    
    val = squeeze(dist(:,:,3,3));
    f = figure;boxplot(val');grid on;set(gca,'ylim',[0 10]);
    f.CurrentAxes.XTickLabel = LandmarkNames(:,1);
    f.CurrentAxes.XTickLabelRotation = 90;
    f.Color = [1 1 1];
    ylabel('(mm)');
    title('AUTOMATIC observer variations');    
end