function out = deleteProperties(Data,index)
         out = Data;
         fullindex = (1:length(Data.Properties.Names));
         index = setdiff(fullindex,index);
         out = rmfield(out,'Properties');
         out.Properties.Names = Data.Properties.Names(index);
         out.Properties.Values = Data.Properties.Values(index,:);
         out.Properties.nr = length(index);
end