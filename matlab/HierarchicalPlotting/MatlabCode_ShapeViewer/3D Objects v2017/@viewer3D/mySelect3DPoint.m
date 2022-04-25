function p = mySelect3DPoint(obj)
         %if nargin<2, slowdown = true;end
         obj.MuteCallbacks = true;
         obj.CurrentShape.PatchHandle.ButtonDownFcn = @shape3D.myClick; 
         %emulateExtentClick(obj);
         %if slowdown, pause(0.05); end
         %pause(0.05);
         obj.MouseRobot.mousePress(obj.MouseExtent);
         pause(0.05);
         p = obj.CurrentShape.PatchClickedPoint';
         obj.CurrentShape.PatchClickedPoint = [];
         obj.CurrentShape.PatchHandle.ButtonDownFcn = [];
         obj.MouseRobot.mouseRelease(obj.MouseExtent);
         %pause(0.05);
         obj.MuteCallbacks = false;
end

% function p = mySelect3DPoint(obj,who)
%          if nargin<2, who = 'currentshape'; end
%          switch who
%              case 'currentshape'
%                 obj.CurrentShape.PatchHandle.ButtonDownFcn = @shape3D.myClick; 
%                 emulateExtentClick(obj);
%                 pause(0.05); 
%                 p = obj.CurrentShape.PatchClickedPoint';
%                 obj.CurrentShape.PatchClickedPoint = [];
%                 obj.CurrentShape.PatchHandle.ButtonDownFcn = [];
%              case 'areaselection'
%                 obj.AreaSelection.PatchHandle.ButtonDownFcn = @shape3D.myClick;
%                 emulateExtentClick(obj);
%                 pause(0.05);
%                 p = obj.AreaSelection.PatchClickedPoint';
%                 obj.AreaSelection.PatchClickedPoint = [];
%                 obj.AreaSelection.PatchHandle.ButtonDownFcn = [];
%              case 'areaselectiondummy'
%                 obj.AreaSelectionDummy.PatchHandle.ButtonDownFcn = @shape3D.myClick;
%                 emulateExtentClick(obj);
%                 pause(0.05);
%                 p = obj.AreaSelectionDummy.PatchClickedPoint';
%                 obj.AreaSelectionDummy.PatchClickedPoint = [];
%                 obj.AreaSelectionDummy.PatchHandle.ButtonDownFcn = [];
%          end
% end