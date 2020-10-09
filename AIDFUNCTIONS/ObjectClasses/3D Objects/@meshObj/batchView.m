function batchView(obj)
         v = viewer(obj,'CalledProcess','Batch Processing');
         obj.Selected = true;
         waitfor(v.Figure);% this is the key point, waiting for the figure to close or press enter
end