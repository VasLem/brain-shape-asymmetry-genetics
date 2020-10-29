function out = getModel(obj,D)
    if nargout == 1,obj = clone(obj);out = obj;end
    if nargin==2&&isstruct(D)
       if isfield(D,'Shape')
            checkShape(obj);
            getModel(obj.Shape,D);
       end
       if isfield(D,'TextureColor')
          checkTexture(obj);
          if strcmp(obj.Texture.Type,'texturePCA'), getModel(obj.Texture,D);end
       end
       if isfield(D,'TextureMap')
          checkTexture(obj,'textureMapPCA');
          if strcmp(obj.Texture.Type,'textureMapPCA'), getModel(obj.Texture,D);end
       end 
    end
    D = obj.WS*eye(obj.nrSC)*(obj.Shape.Tcoeff');
    D = [D;obj.Texture.Tcoeff'];
    getModel@PCA(obj,D);
end