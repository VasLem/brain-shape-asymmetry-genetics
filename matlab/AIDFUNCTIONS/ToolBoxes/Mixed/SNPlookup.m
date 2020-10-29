function [ind12,ind21,raw] = SNPlookup(Array1,Array2,clean)
   % Array1(ind21) = Array2(ind12);
   if nargin < 3, clean = false; end;
   ind12 = nan*ones(1,length(Array1));
   if clean
       Array1 = cleanUpRSArray(Array1);
       Array2 = cleanUpRSArray(Array2);
   end
   [path,ID] = setupParForProgress(length(Array1));
   parfor i=1:length(Array1)
      index = find(strcmp(Array1{i},Array2));
      if ~isempty(index), ind12(i) = index(1); end
      parfor_progress;
   end
   closeParForProgress(path,ID);
   raw = ind12;
   ind21 = find(~isnan(ind12));
   ind12 = ind12(ind21);
end