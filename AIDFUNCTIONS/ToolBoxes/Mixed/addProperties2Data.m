function Data = addProperties2Data(Data,prop,reduce)
         if nargin<3, reduce = false; end
         Data.Properties.Names = prop(1,2:end);
         Data.Properties.nr = length(Data.Properties.Names);
         Data.Properties.Values = cell(Data.Properties.nr,length(Data.Names));
         Included = [];
         for i=1:1:length(Data.Names)
             name = cleanUpString(Data.Names{i});
             for j=2:1:size(prop,1)
                 tmp = cleanUpString(prop{j,1});
                 if ~isempty(strfind(name,tmp))
                     Data.Properties.Values(:,i) = prop(j,2:end);
                     Included = [Included i]; %#ok<AGROW>
                     break;
                 end
             end
         end
         if reduce
             Data = reduceData(Data,Included);
         end
end
