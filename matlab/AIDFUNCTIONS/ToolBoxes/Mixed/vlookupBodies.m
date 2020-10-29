function [ind12,ind2,raw] = vlookupBodies(Array1,Array2)
   ind12 = nan*ones(1,length(Array1));
   parfor i=1:length(Array1)
       name = Array1{i};
       if ~ischar(name), name = num2str(name); end
       name = cleanUpString(name);
       for j=1:1:length(Array2)
           tmp = Array2{j};
           %if strcmp(tmp(1),'M'), tmp
           tmp = tmp(1:3);
           if strcmp(tmp,name)
              ind12(i) = j; 
              break;
           end
       end
   end
   raw = ind12;
   ind2 = find(~isnan(ind12));
   ind12 = ind12(ind2);
end