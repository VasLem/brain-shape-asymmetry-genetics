function [ind12,ind21,raw] = vlookupFast(Array1,Array2,rows)
         % [ind12,ind21,raw] = vlookupFast(Array1,Array2)
         % input = Array1,Array2
         % output = ind12,ind21
         % Array2(ind12) = Array1(ind21);
         if nargin<3, rows = false; end
         if rows
            [~,ind12] = ismember(Array1,Array2,'rows');
         else
            [~,ind12] = ismember(Array1,Array2);
         end
         raw = ind12;
         ind21 = find(ind12);
         ind12 = ind12(ind21);
end