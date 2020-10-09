function D = superimposeLSOntoAverage(obj,D)
         D = getData(obj,D);
         if isempty(obj.Average),getAverage(obj,D);end
         if nargin < 3, normsize = false;iter = 3;end
         if nargin < 4, iter = 3; end
         if normsize 
            T = scaledRigidTM;
         else
            T = rigidTM;
         end
         nrSamples = size(D,2);
         nrV = obj.Average.nrV;
         f1 = waitbar(0,'LSGP');drawnow
         drawnow;
         for j=1:1:iter
            %centerVertices(obj.Average);
            alignWithAxes(obj.Average);
            avgpoints = obj.Average.Vertices;
            if normsize, avgpoints = avgpoints/centroidSize(obj.Average); end   
            parfor i=1:nrSamples
                Tfor = clone(T);
                points = reshape(D(:,i),3,nrV);
                if ~normsize,Tfor.c = [];end
                match(Tfor,avgpoints,points);
                eval(Tfor,points); %#ok<PFBFN>
                D(:,i) = Tfor.Evaluation.Vertices(:);
                delete(Tfor);
            end
            waitbar((j/iter),f1);drawnow;
            drawnow;
            getAverage(obj,D);
         end
         delete(T);
         delete(f1);drawnow;
end