function out = getProcrustesDistance(S1,S2,W)
         if nargin < 3, W = []; end
         if size(S1,1)==1||size(S1,2)==1, S1 = reshape(S1(:),3,length(S1)/3);end
         if size(S2,1)==1||size(S2,2)==1, S2 = reshape(S2(:),3,length(S2)/3);end
         SS = (S2-S1).^2;
         Distances = sqrt(sum(SS));
         RMSE = sqrt(mean(Distances.^2));
         PD = sum(SS(:));
         out.SS = SS;
         out.Distances = Distances;
         out.RMSE = RMSE;
         out.PD = PD;
         if isempty(W), return; end
         FW = repmat(W,3,1);
         WSS = FW.*SS;
         WPD = sum(WSS(:));
         WRMSE = sqrt(sum(W.*(Distances).^2)/sum(W));
         out.WPD = WPD;
         out.WRMSE = WRMSE;
end