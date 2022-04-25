function createPatch(obj)
         if isempty(obj.Axes)||~ishandle(obj.Axes), obj.Axes = gca; end
         obj.ph = patch('Parent',obj.Axes,'Visible','off','UserData',obj);
         setPatch(obj,'All');
         updateChildren(obj,'Patch Created');
end % create Patch handle.