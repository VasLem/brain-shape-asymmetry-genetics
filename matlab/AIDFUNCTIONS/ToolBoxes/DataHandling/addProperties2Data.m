function [Data,Missed] = addProperties2Data(Data,prop,reduce)
         if nargin<3, reduce = false; end
         Data.Properties.Names = prop(1,2:end);
         Data.Properties.nr = length(Data.Properties.Names);
         Data.Properties.Values = cell(Data.Properties.nr,length(Data.Names));
         Data.Properties.Values(:) = {nan};
         Included = [];
         IncProp = ones(1,size(prop,1));
         IncProp(1) = 0;
         for i=1:1:length(Data.Names)
             name = cleanUpString(Data.Names{i});
             for j=2:1:size(prop,1)
                 tmp = prop{j,1};
                 if ~ischar(tmp), tmp = num2str(tmp); end
                 tmp = cleanUpString(tmp);
                 if ~isempty(strfind(name,tmp))
                     Data.Properties.Values(:,i) = prop(j,2:end);
                     Included = [Included i]; %#ok<AGROW>
                     IncProp(j) = 0;
                     break;
                 end
             end
         end
         if reduce
            Data = reduceData(Data,Included);
         end
         Missed = prop(find(IncProp),1);    
end
