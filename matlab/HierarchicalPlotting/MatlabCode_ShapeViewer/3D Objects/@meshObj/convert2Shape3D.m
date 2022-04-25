function out = convert2Shape3D(obj)
    out = shape3D;
    out.Vertices = obj.Vertices';
    out.Faces = obj.Faces';
    out.VertexValue = obj.Value;
    if ~isempty(obj.UV), out.VertexRGB = obj.TextureColor;end
    out.UV = obj.UV;
    if ~isempty(obj.TextureMap), out.TextureMap = obj.TextureMap.Image; end
    out.SingleColor = obj.SingleColor;
    out.Material = obj.Material;
    out.ViewMode = obj.ViewMode;
    out.ColorMode = obj.ColorMode;
    out.VertexSize = obj.MarkerSize;
end