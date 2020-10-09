function D = LSGenProcrustes(obj,D,normsize,iter,refscan)
         D = getData(obj,D);
         if isempty(obj.Average),getAverage(obj,D);end
         if nargin < 3, normsize = false;iter = 3;end
         if nargin < 4, iter = 3; end
         if nargin < 5, refscan = []; end
         if normsize 
            T = scaledRigidTM;
         else
            T = rigidTM;
         end
         nrSamples = size(D,2);
         nrV = obj.Average.nrV;
         %f1 = waitbar(0,'LSGP');drawnow
         %drawnow;
         for j=1:1:iter
            disp(['ITER ' num2str(j) ' / ' num2str(iter)]);
            disp('');
            if ~isempty(refscan) 
               alignPose(obj.Average,refscan,0);
            else
               %centerVertices(obj.Average);
               alignWithAxes(obj.Average);
            end
            avgpoints = obj.Average.Vertices;
            if normsize, avgpoints = avgpoints/centroidSize(obj.Average); end
            [path,ID] = setupParForProgress(nrSamples);
            parfor i=1:nrSamples
                Tfor = clone(T);
                points = reshape(D(:,i),3,nrV);
                if ~normsize,Tfor.c = [];end
                match(Tfor,avgpoints,points);
                eval(Tfor,points); %#ok<PFBFN>
                D(:,i) = Tfor.Evaluation.Vertices(:);
                parfor_progress;
            end
            closeParForProgress(path,ID);
            getAverage(obj,D);  
         end
end