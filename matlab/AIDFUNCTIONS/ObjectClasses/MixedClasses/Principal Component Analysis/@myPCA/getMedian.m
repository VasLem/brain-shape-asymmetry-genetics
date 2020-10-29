function out = getMedian(obj,in)
% out = getMedian(obj,in) or getAverage(obj,in)
% Compute the element-wise Median entity based on give training data
% INPUT
% obj = PCA space object
% in = Data matrix or Data structure
% OUTPUT
% out = element-wise median entity
%
% created by Peter Claes
         in = getData(obj,in);
         med = nanmedian(in,2);
         out = Vec2Struc(obj,med);
end