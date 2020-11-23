function [AlignedLandmarks,AvgLandmarks,CentroidSizes,v] = GeneralizedProcrustesAnalysis(Landmarks,TemplateLandmarks,iter,scale,reflect,progress,display)
         if nargin<7, display = false;end
         if nargin<6, progress = false;end
         if nargin<5, reflect = true;end
         if nargin<4, scale = true;end
         if nargin<3, iter = 3; end
         if ismatrix(Landmarks)% turn matrix into 3 dim array
            N = size(Landmarks,2);
            nLM = size(Landmarks,1)/3;
            Landmarks = reshape(Landmarks,nLM,3,N);
         else
            N = size(Landmarks,3);
            nLM = size(Landmarks,1);
         end
         if ~isempty(TemplateLandmarks)
            AvgLandmarks = TemplateLandmarks.Vertices;
         else
            AvgLandmarks = mean(Landmarks,3);
         end
         AlignedLandmarks = Landmarks;
         avg = mean(AlignedLandmarks,1);
        Differences = avg - AlignedLandmarks;
        Distances = sqrt(sum(Differences.^2,2));
        CentroidSizes = sqrt(sum(Distances.^2)/nLM);      
         for i=1:iter
             if progress, disp([num2str(i) ' out of ' num2str(iter)]);end
%              if progress, [path,ID] = setupParForProgress(N);end
%              parfor n=1:1:N
              for n=1:1:N 
                tmp = squeeze(AlignedLandmarks(:,:,n));
                 [~,AlignedLandmarks(:,:,n),~] = procrustes(AvgLandmarks,tmp,'Scaling',scale,'Reflection',reflect);
%                  if progress, parfor_progress;end
             end
%              if progress, closeParForProgress(path,ID);end
             AvgLandmarks = mean(AlignedLandmarks,3);
             if scale
                tmp = shape3D;tmp.Vertices = AvgLandmarks;
                AvgLandmarks = AvgLandmarks./sqrt(sum(AvgLandmarks.^2,'all'));
             end
         end
         v = [];
         if display
            display3DLandmarks(TemplateLandmarks, AvgLandmarks);
         end
         if ~isempty(TemplateLandmarks)
             tmp = AvgLandmarks;
             AvgLandmarks = clone(TemplateLandmarks);
             AvgLandmarks.Vertices = tmp;
         end
end