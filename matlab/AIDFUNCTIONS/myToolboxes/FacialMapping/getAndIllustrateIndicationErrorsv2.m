function [intra,inter] = getAndIllustrateIndicationErrorsv2(DATA,Template,AvgLM,LandmarkNames)
    
     intra = [];
     inter = [];

     N = length(DATA);
     % intra observer indication errors
     Template = clone(Template);
     AvgLM = clone(AvgLM);
     v = viewer(Template);
     viewer(AvgLM,v); 
     dist = zeros(AvgLM.nVertices,N,3,3,2);
     Vertices = [];
     %i=1;o=1;t=1;
     for i=1:1:N
        for o=1:2
            for t=1:3
                for tt=1:3
                    diff = DATA{i}.LM{o,t}.Vertices - DATA{i}.LM{o,tt}.Vertices;
                    dist(:,i,t,tt,o) = sqrt(sum(diff.^2,2));
                    Vertices = [Vertices; AvgLM.Vertices+diff]; %#ok<AGROW>
                end
            end
        end
     end
     IntraLM = clone(AvgLM);
     IntraLM.Vertices = Vertices;
     IntraLM.SingleColor = [1 0 0];
     IntraLM.ColorMode = 'single';
     IntraLM.VertexSize = 10;
     viewer(IntraLM,v);
     
     nLM = AvgLM.nVertices;
     DIST = zeros(41*2*3*3,nLM);
     for l=1:1:nLM
         tmp = squeeze(dist(l,:,:,:));
         DIST(:,l) = tmp(:);
     end
     index = find(DIST(:,1));
     DIST = DIST(index,:);
     f = figure;boxplot(DIST);grid on;set(gca,'ylim',[0 10]);
     f.CurrentAxes.XTickLabel = LandmarkNames(:,1);
     f.CurrentAxes.XTickLabelRotation = 90;
     f.Color = [1 1 1];
     ylabel('(mm)');
     title('INTRA observer variations');
     
     Template = clone(Template);
     AvgLM = clone(AvgLM);
     v = viewer(Template);
     viewer(AvgLM,v); 
     dist = zeros(AvgLM.nVertices,N,3,3,2);
     Vertices = [];
     %i=1;o=1;t=1;
     for i=1:1:N
        for o=1:2
            oo = setdiff(1:2,o);
            for t=1:3
                for tt=1:3
                    diff = DATA{i}.LM{o,t}.Vertices - DATA{i}.LM{oo,tt}.Vertices;
                    dist(:,i,t,tt,o) = sqrt(sum(diff.^2,2));
                    Vertices = [Vertices; AvgLM.Vertices+diff]; %#ok<AGROW>
                end
            end
        end
     end
     InterLM = clone(AvgLM);
     InterLM.Vertices = Vertices;
     InterLM.SingleColor = [1 0 0];
     InterLM.ColorMode = 'single';
     InterLM.VertexSize = 10;
     viewer(InterLM,v);
     
     nLM = AvgLM.nVertices;
     DIST = zeros(41*2*3*3,nLM);
     for l=1:1:nLM
         tmp = squeeze(dist(l,:,:,:));
         DIST(:,l) = tmp(:);
     end
     index = find(DIST(:,1));
     DIST = DIST(index,:);
     f = figure;boxplot(DIST);grid on;set(gca,'ylim',[0 10]);
     f.CurrentAxes.XTickLabel = LandmarkNames(:,1);
     f.CurrentAxes.XTickLabelRotation = 90;
     f.Color = [1 1 1];
     ylabel('(mm)');
     title('INTER observer variations');
end