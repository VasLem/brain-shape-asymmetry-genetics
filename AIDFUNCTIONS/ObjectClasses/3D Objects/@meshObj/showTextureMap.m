function showTextureMap(obj,varargin)
    figure;imshow(obj.TextureMap);hold on;
    if isempty(find(strcmp(varargin,'UV'))), return; end
    uv = obj.UV;
    uv(1,:) = uv(1,:)*obj.TextureMap.ImageSize(1);
    uv(2,:) = uv(2,:)*obj.TextureMap.ImageSize(2);
    fastrbf_view(uv);
end
