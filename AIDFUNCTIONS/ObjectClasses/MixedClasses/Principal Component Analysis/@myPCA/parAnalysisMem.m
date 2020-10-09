function out = parAnalysisMem(obj,D,alpha,nShuffle,display)
% Parallel Analysis to determine the amount of principle components needed
%
%
% created by Peter Claes
% compute the original model
         obj = clone(obj);
         getAverage(obj,D);
         D = D-repmat(obj.AvgVec,1,size(D,2));
         obj.Centering = false;
         getModel(obj,D);
         latent = obj.EigVal;
         if nShuffle==0, out = []; return; end
         latentShuffle = zeros(length(latent), nShuffle);
         x = D';
         [path,ID] = setupParForProgress(nShuffle);
         for iShuffle = 1:nShuffle
                 rng(iShuffle);
                 xShuffle = x;
                 for dim = 1:size(x,2)
                         xShuffle(:, dim) = x(randperm(size(x,1)), dim);
                 end
                 forobj = clone(obj);
                 %getAverage(forobj,xShuffle');
                 getModel(forobj,xShuffle');
                 latentShuffle(:,iShuffle) = forobj.EigVal;
                 parfor_progress;
         end
         closeParForProgress(path,ID);
         latentHigh = quantile(latentShuffle', 1-alpha);
         latentLow = quantile(latentShuffle', alpha);
         latentMean = mean(latentShuffle,2);
         latentMedian = median(latentShuffle,2);
         out.latent = latent;
         out.latentHigh = latentHigh;
         out.latentLow = latentLow;
         out.latentMean = latentMean;
         out.nrPClow = sum(latent>=latentLow');
         out.nrPChigh = sum(latent>=latentHigh');
         out.nrPCmean = sum(latent>=latentMean);
         out.nrPCmedian = sum(latent>=latentMedian);
         if display
                plot(latent,'bo-');
                hold on;
                plot(latentHigh,'ro-');
                plot(latentLow,'go-');
                plot(latentMean,'k--');
                hold off;
                legend('data',sprintf('shuffled %d%%', 100*(1-alpha)),sprintf('shuffled %d%%', 100*alpha));
         end
end

