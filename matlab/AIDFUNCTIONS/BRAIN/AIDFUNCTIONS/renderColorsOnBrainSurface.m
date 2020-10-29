function scan = renderColorsOnBrainSurface(rend,colors,peer)
    scan = clone(rend.RefScan);
    scan.VertexRGB = colors;
    scan.ColorMode = 'Texture';
    scan.Material = 'Dull';
    scan.ViewMode = 'solid';
    scan.RenderAxes = peer;
    scan.Visible = true;
    scan.PatchHandle.FaceColor = 'flat';
    axis(peer,'image');
    axis(peer,'off');
    %colorbar(peer,'NorthOutside');
    colorbar(peer,'off');   
end