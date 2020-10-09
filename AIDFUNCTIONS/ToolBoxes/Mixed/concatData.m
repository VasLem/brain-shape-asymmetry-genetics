function out = concatData(Data,Data2)
         if isfield(Data,'Names'),out.Names = [Data.Names, Data2.Names];end
         if isfield(Data,'Shape'),out.Shape = [Data.Shape, Data2.Shape];end
         if isfield(Data,'TextureColor'),out.TextureColor = [Data.TextureColor, Data2.TextureColor];end
         if isfield(Data,'TextureMap'),out.TextureMap = [Data.TextureMap, Data2.TextureMap];end
         if isfield(Data,'Properties'),
             out.Properties.Names = Data.Properties.Names;
             out.Properties.nr = Data.Properties.nr;
             out.Properties.Values = [Data.Properties.Values, Data2.Properties.Values];
         end
end