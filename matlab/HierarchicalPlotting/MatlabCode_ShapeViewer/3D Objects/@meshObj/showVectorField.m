function v = showVectorField(startscan,endscan,v)
         if nargin < 3
            v = viewer3DObj; 
         end
         VF = vectorField3D;
         VF.StartPoints = startscan.Vertices;
         VF.EndPoints = endscan.Vertices;
         viewer(VF,'Viewer',v);
         shell = clone(startscan);
         shell.Axes = v.RenderAxes;
         shell.SingleColor = [0.5 0.5 0.5];
         if ~isempty(shell.PoseLM), delete(shell.PoseLM); end
%          reduceTriangles(shell,0.2);
         shell.ViewMode = 'Wireframe';
         shell.MarkerSize = 10;
         shell.ColorMode = 'Single';
         shell.Visible = true;
         if ~isempty(shell.Border), shell.Border.Visible = false; end
         v.Renderer = 'opengl';
         shell.Alpha = 0.6;
         v.SceneLightVisible = true;
         v.SceneLightLinked = true;
         shell.Material = 'Dull';
         set(VF.ph,'LineWidth',2);
end