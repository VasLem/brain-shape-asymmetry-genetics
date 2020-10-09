function [intra,inter] = getAndIllustrateIndicationErrors(DATA,Template,AvgLM,LandmarkNames)


     N = length(DATA);
     % intra observer indication errors
     Template = clone(Template);
     AvgLM = clone(AvgLM);
     v = viewer(Template);
     viewer(AvgLM,v);

     
     dist = zeros(AvgLM.nVertices,N,2,3);
     Vertices = [];
     %i=1;o=1;t=1;
     for i=1:1:N
        for o=1:2
            for t=1:3
                diff = DATA{i}.LM{o,t}.Vertices - DATA{i}.AvgLM{o}.Vertices;
                dist(:,i,o,t) = sqrt(sum(diff.^2,2));
                Vertices = [Vertices; AvgLM.Vertices+diff]; %#ok<AGROW>
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
     DIST = zeros(41*2*3,nLM);
     for l=1:1:nLM
         tmp = squeeze(dist(l,:,:,:));
         DIST(:,l) = tmp(:);
     end
     f = figure;boxplot(DIST);grid on;set(gca,'ylim',[0 10]);
     f.CurrentAxes.XTickLabel = LandmarkNames(:,1);
     f.CurrentAxes.XTickLabelRotation = 90;
     f.Color = [1 1 1];
     ylabel('(mm)');
     title('INTRA observer variations');
     
     intra.MeanDistances = mean(squeeze(mean(dist,2)),3);
     intra.LM = clone(IntraLM);
     intra.AvgLM = AvgLM;
     intra.Distances = dist;
     
     disp(intra.MeanDistances);
   
     
%      % inter observer indiciation errors
%      Template = clone(Template);
%      AvgLM = clone(AvgLM);
%      v = viewer(Template);
%      viewer(AvgLM,v);
%       
%      
%      dist = zeros(AvgLM.nVertices,N,2,3);
%      Vertices = [];
%      %i=1;o=1;t=1;
%      for i=1:1:N
%         for o=1:2
%             oo = setdiff(1:2,o);
%             for t=1:3
%                 diff = DATA{i}.LM{o,t}.Vertices - DATA{i}.AvgLM{oo}.Vertices;
%                 dist(:,i,o,t) = sqrt(sum(diff.^2,2));
%                 Vertices = [Vertices; AvgLM.Vertices+diff]; %#ok<AGROW>
%             end
%         end
%      end
%      InterLM = clone(AvgLM);
%      InterLM.Vertices = Vertices;
%      InterLM.SingleColor = [0 0 1];
%      InterLM.ColorMode = 'single';
%      InterLM.VertexSize = 10;
%      viewer(InterLM,v);
%      
%      inter.MeanDistances = mean(squeeze(mean(dist,2)),3);
%      inter.LM = clone(InterLM);
%      inter.AvgLM = AvgLM;
%      inter.Distances = dist;
%      
%      disp(inter.MeanDistances);

     % inter observer indiciation errors
     Template = clone(Template);
     AvgLM = clone(AvgLM);
     v = viewer(Template);
     viewer(AvgLM,v);
      
     
     dist = zeros(AvgLM.nVertices,N,2);
     Vertices = [];
     %i=1;o=1;t=1;
     for i=1:1:N
        for o=1:2
            oo = setdiff(1:2,o);
            diff = DATA{i}.AvgLM{o}.Vertices - DATA{i}.AvgLM{oo}.Vertices;
            Vertices = [Vertices; AvgLM.Vertices+diff]; %#ok<AGROW>
            dist(:,i,o) = sqrt(sum(diff.^2,2));
%             for t=1:3
%                 diff = DATA{i}.LM{o,t}.Vertices - DATA{i}.AvgLM{oo}.Vertices;
%                 dist(:,i,o,t) = sqrt(sum(diff.^2,2));
%                 
%             end
        end
     end
     InterLM = clone(AvgLM);
     InterLM.Vertices = Vertices;
     InterLM.SingleColor = [0 0 1];
     InterLM.ColorMode = 'single';
     InterLM.VertexSize = 10;
     viewer(InterLM,v);
     
     inter.MeanDistances = mean(squeeze(mean(dist,2)),2);
     inter.LM = clone(InterLM);
     inter.AvgLM = AvgLM;
     inter.Distances = dist;
     
     disp(inter.MeanDistances);


     % inter observer indiciation errors
     Template = clone(Template);
     AvgLM = clone(AvgLM);
     v = viewer(Template);
     viewer(AvgLM,v);
     viewer(intra.LM,v);
     viewer(inter.LM,v);
     
     DIST = zeros(41*2,nLM);
     for l=1:1:nLM
         tmp = squeeze(dist(l,:,:));
         DIST(:,l) = tmp(:);
     end
     f = figure;boxplot(DIST);grid on;set(gca,'ylim',[0 10]);
     f.CurrentAxes.XTickLabel = LandmarkNames(:,1);
     f.CurrentAxes.XTickLabelRotation = 90;
     f.Color = [1 1 1];
     ylabel('(mm)');
     title('INTER observer variations');





end