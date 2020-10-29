function [DL,DR,D] = LSGenAsymProcrustes(obj,DL,DR,normsize,iter)
         DL = getData(obj,DL);
         DR = getData(obj,DR);
         D = (DL+DR)/2;
         if isempty(obj.Average),getAverage(obj,D);end
         if nargin < 3, normsize = false;iter = 3;end
         if nargin < 4, iter = 3; end
         if normsize 
            T = scaledRigidTM;
         else
            T = rigidTM;
         end
         nrSamples = size(D,2);
         f1 = statusbar('LSAsymGP');
         drawnow;
         for j=1:1:iter
            centerVertices(obj.Average);
            avgpoints = obj.Average.Vertices;
            if normsize, avgpoints = avgpoints/centroidSize(obj.Average); end
            nrV = obj.Average.nrV;
            parfor i=1:nrSamples
                Tfor = clone(T);
                pointsL = reshape(DL(:,i),3,nrV);
                pointsR = reshape(DR(:,i),3,nrV);
                pointsLR = (pointsL+pointsR)/2;
                if ~normsize,Tfor.c = [];end
                match(Tfor,avgpoints,pointsLR);
                eval(Tfor,pointsL);
                DL(:,i) = Tfor.Evaluation.Vertices(:);
                eval(Tfor,pointsR);
                DR(:,i) = Tfor.Evaluation.Vertices(:);
                delete(Tfor);
            end
            statusbar((j/iter),f1);
            drawnow;
            D = (DL+DR)/2;
            getAverage(obj,D);
         end
         delete(T);
         delete(f1);
end