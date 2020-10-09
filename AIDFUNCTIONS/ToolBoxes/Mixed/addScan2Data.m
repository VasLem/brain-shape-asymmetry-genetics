function out = addScan2Data(Data,scan)
         out = Data;
         if isfield(Data,'Names'),out.Names{end+1} = scan.Tag;end
         if isfield(Data,'Shape'),out.Shape(:,end+1) = scan.Vertices(:);end
         if isfield(Data,'TextureColor'),out.TextureColor(:,end+1) = scan.TextureColor(:);end
         if isfield(Data,'TextureMap'),out.TextureMap(:,end+1) = scan.TextureMap.uint8Image(:);end
         if isfield(Data,'Properties')
             out.Properties.Values(:,end+1) = [scan.Prop(:); scan.OrthoProp.TriangularMeasures(:); scan.OrthoProp.StableLM(:)];
         end
end