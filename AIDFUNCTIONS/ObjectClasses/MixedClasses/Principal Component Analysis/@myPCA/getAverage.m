function out = getAverage(obj,in)
% out = getAverage(obj,in) or getAverage(obj,in)
% Compute the Average for the PCA space based on give training data
% INPUT
% obj = PCA space object
% in = Data matrix or Data structure
% OUTPUT
% out = new PCA space object or out = updated obj;
%
% created by Peter Claes   
         in = getData(obj,in);
         if obj.Centering
            gem = nanmean(in,2);
         else
            gem = zeros(size(in,1),1);
         end
         if nargout == 1, out = gem; return; end
         obj.AvgVec = gem;
end