function D = analyseShapeAlignment(obj,D)
         
         if isempty(obj.Average),getAverage(obj,D);end
         nrSamples = size(D,2);     
         iter = 1;
         i = 150;
         %v = viewer(obj.Average);
         face = clone(obj.Average);
         face.SingleColor = [0.5 0.5 0.5];
         face.Vertices = reshape(D(:,i),3,obj.Average.nrV);
         face.Value = vDistances(obj.Average,face);
         face.ColorMode = 'Indexed';
         viewer(face);
         
%          
%          for j=1:1:iter
%             avgpoints = obj.Average.Vertices;
%             f1 = statusbar(['Alignment ' num2str(j) '/' num2str(iter)]);
%             for i=1:1:nrSamples
%                 points = reshape(D(:,i),3,obj.Average.nrV);
%                 if ~scale,T.c = [];end
%                 match(T,avgpoints,points);
%                 eval(T,points);
%                 D(:,i) = T.Evaluation.Vertices(:);
%                 clear points;
%                 statusbar((i/nrSamples),f1);
%             end
%             delete(f1);
%             getAverage(obj,D);
%          end
%          delete(T);
end