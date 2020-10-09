function [ind12,ind21,raw] = vlookup(Array1,Array2,clean)
   % Array1(ind21) = Array2(ind12);
   if nargin < 3, clean = true; end;
   ind12 = nan*ones(1,length(Array1));
   if clean
       parfor i=1:length(Array1)
           name = Array1{i};
           if ~ischar(name), name = num2str(name); end
           Array1{i} = cleanUpString(name);
       end
       parfor i=1:length(Array2)
           name = Array2{i};
           if ~ischar(name), name = num2str(name); end
           Array2{i} = cleanUpString(name);
       end
   end
   [path,ID] = setupParForProgress(length(Array1));
   parfor i=1:length(Array1)
      tmp = strcmpi(Array1{i},Array2);
      index = find(tmp);
      if ~isempty(index), ind12(i) = index(1); end
      parfor_progress;
   end
   closeParForProgress(path,ID);
   raw = ind12;
   ind21 = find(~isnan(ind12));
   ind12 = ind12(ind21);
end