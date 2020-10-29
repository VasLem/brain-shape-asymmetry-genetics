function renderPatch(obj)
         obj.Visible =  true;
         if isempty(obj.ph)||~ishandle(obj.ph), createPatch(obj); return; end
end % render the object.