function getAverage(obj,D)
         if ~isstruct(D), return; end
         if isfield(D,'Shape')
            checkShape(obj);
            getAverage(obj.Shape,D);
         end
         if isfield(D,'TextureColor')
            checkTexture(obj);
            if strcmp(obj.Texture.Type,'texturePCA'), getAverage(obj.Texture,D);end
         end
         if isfield(D,'TextureMap')
            checkTexture(obj,'textureMapPCA');
            if strcmp(obj.Texture.Type,'textureMapPCA'), getAverage(obj.Texture,D);end
         end
end