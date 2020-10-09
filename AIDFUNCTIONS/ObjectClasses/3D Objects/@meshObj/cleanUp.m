function cleanUp(obj)
         %oldVis = obj.Visible;
         %oldAxes = obj.Axes;
         v = viewer(obj,'CalledProcess','Batch Processing');
         obj.Selected = true;
         waitfor(v.Figure);     
end