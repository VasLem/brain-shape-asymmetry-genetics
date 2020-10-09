function D = alignShape(obj,D,scale,index)
         D = getData(obj,D);
         if isempty(obj.Average),getAverage(obj,D);end
         if nargin < 3, scale = false;index = [];end
         if nargin < 4, index = []; end
         if scale
            T = scaledRigidTM;
         else
            T = rigidTM;
         end
         if isempty(index), index = (1:obj.Average.nrV); end
         nrSamples = size(D,2);     
         iter = 1;         
         for j=1:1:iter
            avgpoints = obj.Average.Vertices;
            f1 = statusbar(['Alignment ' num2str(j) '/' num2str(iter)]);
            for i=1:1:nrSamples
                points = reshape(D(:,i),3,obj.Average.nrV);
                if ~scale,T.c = [];end
                match(T,avgpoints(:,index),points(:,index));
                eval(T,points);
                D(:,i) = T.Evaluation.Vertices(:);
                clear points;
                statusbar((i/nrSamples),f1);
                drawnow;
            end
            delete(f1);
            getAverage(obj,D);
         end
         delete(T);
end