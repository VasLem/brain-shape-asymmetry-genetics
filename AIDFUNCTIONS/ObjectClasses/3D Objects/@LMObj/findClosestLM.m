function [index] = findClosestLM(obj,p,min_dist)
         if isempty(obj.Location), index = []; return; end
         [distances,I] = sort(sqrt(sum((repmat(p,1,size(obj.Location,2))-obj.Location).^2)));
         if distances(1) <= min_dist
            index = I(1);
         else
            index = [];
         end
end
       