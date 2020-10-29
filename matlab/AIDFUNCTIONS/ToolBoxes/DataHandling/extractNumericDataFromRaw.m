function out = extractNumericDataFromRaw(in)
         out = in;
         for i=1:1:numel(out)
             tmp = out{i};
             if ~isnumeric(tmp)
                 tmp = nan;
             end
             out{i} = tmp;
         end
         out = cell2mat(out);
end