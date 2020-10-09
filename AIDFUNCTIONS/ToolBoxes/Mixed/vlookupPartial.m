function [ind12,ind21,raw] = vlookupPartial(Array1,Array2)
   ind12 = nan*ones(1,length(Array1));
   [path,ID] = setupParForProgress(length(Array1));
   parfor i=1:length(Array1)
      tmp = strfind(Array2,Array1{i});
      index = 0;
      for j=1:1:length(Array2);
          if isempty(tmp{j}), continue;end
          index = j;
          break;
      end
      if ~(index==0), ind12(i) = index; end
      parfor_progress;
   end
   closeParForProgress(path,ID);
   raw = ind12;
   ind21 = find(~isnan(ind12));
   ind12 = ind12(ind21);
end


% tmp = strfind(RS,testrs);
