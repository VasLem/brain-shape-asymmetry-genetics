function [map,convertuv] = convert2SingleMap(maps)
         nMaps = length(maps);
         good = ones(1,nMaps);
         for i=1:1:nMaps
             if isempty(maps{i}), good(i) = 0; end        
         end
         
         for i=1:1:good
         
         



end