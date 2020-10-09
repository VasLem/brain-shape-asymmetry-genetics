function out = getNumericProperties(prop)
         names = [];
         values = [];
         for r=1:1:prop.nr
             disp(num2str(r));
             tmp = prop.Values(r,:);
             if ischar(tmp{1})% we are dealing with string properties and need to convert them
                un = {};
                for i=1:1:length(tmp)
                       try
                         un = unique([un(:) ; tmp(i)]);
                       catch
                       end
                end
                for i=1:1:length(un)
                    names = [names {[prop.Names{r} ' ' un{i}]}]; %#ok<AGROW>
                    values = [values; strcmp(un{i},tmp)];  %#ok<AGROW>
                end
             else
                 %try
                   values = [values;cell2mat(tmp)]; %#ok<AGROW>
                 %catch
                 %end
                names = [names prop.Names(r)]; %#ok<AGROW>
             end                    
         end
         out.Names = names;
         out.Values = values;
         out.nr = size(values,1);
end