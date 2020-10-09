function indicatePoseLM(obj)
         cobj = clone(obj);
         cobj.Selected = true;
         v = viewer(cobj,'CalledProcess','indicatePoseLM');
         %%
% %          v.CameraPosition = [-563.5361 -30.1358 77.4830];
% %          v.CameraTarget = [3.0005 1.0124 -0.4743];
% %          v.CameraUpVector = [0.1314 0.0829 0.9879];
% %          v.CameraViewAngle = 6.0964;
%          v.DataAspectRatio = [1 1 1];
%          v.SceneLightMode = 'infinite';
%          v.SceneLightPosition = [-563.5361 -30.1358 77.4830];
         %%
         set(v.Toolbar.light_toggle,'State','on');
         set(v.Toolbar.link_toggle,'State','on');
         cobj.Selected = true;
         waitfor(v.Figure);
         %if size(cobj.PoseLM.Vertices,2)==5
             obj.PoseLM = LMObj('Vertices',cobj.PoseLM.Vertices,'Axes',obj.Axes);
             obj.PoseLM.Visible = obj.Visible;
         %end
         delete(cobj);
end