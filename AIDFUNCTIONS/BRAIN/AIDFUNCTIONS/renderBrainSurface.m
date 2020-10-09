function scan = renderBrainSurface(rend,values,peer)
    scan = clone(rend.RefScan);
    scan.VertexValue = values;
    scan.ColorMode = "Indexed";
    scan.Material = 'Dull';
    scan.ViewMode = 'solid';
    scan.RenderAxes = peer;
    scan.Visible = true;
    scan.PatchHandle.FaceColor = 'flat';
    axis(peer,'image');
    axis(peer,'off');
    colorbar(peer,'off');   
end