function out = getNumericProperties(prop)
         names = [];
         values = [];
         for r=1:1:prop.nr
             tmp = prop.Values(r,:);
             if ischar(tmp{1})% we are dealing with string properties and need to convert them
                un = unique(tmp);
                for i=1:1:length(un)
                    names = [names un(i)]; %#ok<AGROW>
                    values = [values; strcmp(un{i},tmp)];  %#ok<AGROW>
                end
             else
                values = [values;cell2mat(tmp)]; %#ok<AGROW>
                names = [names prop.Names(r)]; %#ok<AGROW>
             end                    
         end
         out.Names = names;
         out.Values = values;
         out.nr = size(values,1);
end