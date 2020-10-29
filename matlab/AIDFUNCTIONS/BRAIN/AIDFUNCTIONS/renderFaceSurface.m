function scan = renderFaceSurface(rend,values,peer)
    scan = clone(rend.RefScan);
    scan.VertexValue = values;
    scan.ColorMode = "Indexed";
    scan.Material = 'Dull';
    scan.ViewMode = 'solid';
    scan.RenderAxes = peer;
    
    theta = -pi/4;
    ROT = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    scan.Vertices = (ROT*scan.Vertices')';

    scan.Visible = true;
    %scan.PatchHandle.FaceColor = 'flat';
    axis(peer,'image');
    axis(peer,'off');
    %colorbar(peer,'NorthOutside');
    colorbar(peer,'off');   
end