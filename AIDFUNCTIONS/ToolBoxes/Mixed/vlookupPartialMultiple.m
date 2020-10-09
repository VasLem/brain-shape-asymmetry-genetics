function [ind12] = vlookupPartialMultiple(Array1,Array2)
   ind12 = cell(1,length(Array1));
   [path,ID] = setupParForProgress(length(Array1));
   parfor i=1:length(Array1)
      tmp = strfind(Array2,Array1{i});
      index = [];
      for j=1:1:length(Array2);
          if isempty(tmp{j}), continue;end
          index = [index j];
          break;
      end
      if ~isempty(index), ind12{i} = index; end
      parfor_progress;
   end
   closeParForProgress(path,ID);
%    raw = ind12;
%    ind21 = find(~isnan(ind12));
%    ind12 = ind12(ind21);
end


% tmp = strfind(RS,testrs);
