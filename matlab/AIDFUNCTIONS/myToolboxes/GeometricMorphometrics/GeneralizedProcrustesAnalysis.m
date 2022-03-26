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
         
         avg = mean(Landmarks,1);
         CentroidSizes = sqrt(sum(sqrt(sum((avg - Landmarks).^2,2)).^2)/nLM);
         AlignedLandmarks = mat2cell(Landmarks, size(Landmarks,1), size(Landmarks,2), ones(1,size(Landmarks,3)));
%         Differences = avg - AlignedLandmarks;
%         Distances = sqrt(sum(Differences.^2,2));
        
        
        
         for i=1:iter
             if progress, disp([num2str(i) ' out of ' num2str(iter)]);end
             if progress, ppb=ParforProgressbar(N);end
             parfor n=1:N
%               for n=1:1:N 
                 [~,AlignedLandmarks{n},~] = procrustes(AvgLandmarks,squeeze(AlignedLandmarks{n}),'Scaling',scale,'Reflection',reflect);
                 if progress, ppb.increment();end
             end
             if progress, delete(ppb); end
             AvgLandmarks = AlignedLandmarks{1};
             for n=1:N
                AvgLandmarks = AvgLandmarks + AlignedLandmarks{n};
             end
             AvgLandmarks = AvgLandmarks/N;
             if scale
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
         AlignedLandmarks = cell2mat(AlignedLandmarks);
end