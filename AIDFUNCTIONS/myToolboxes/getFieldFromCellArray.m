function out = getFieldFromCellArray(in,fieldnames)
                out = [];
                if isempty(in), return; end
                nrCells = length(in);
                nFields = length(fieldnames);
                for f=1:1:nFields
                    eval(['out.' fieldnames{f} ' = [];']);
                end
                for i=1:1:nrCells
                   if isempty(in{i})
                      for f=1:1:nFields
                        eval(['out.' fieldnames{f} ' = [out.' fieldnames{f} ' nan];']);
                      end
                   else
                      for f=1:1:nFields
                        eval(['out.' fieldnames{f} ' = [out.' fieldnames{f} ' in{i}.' fieldnames{f} '];']);
                      end
                   end
                end
end