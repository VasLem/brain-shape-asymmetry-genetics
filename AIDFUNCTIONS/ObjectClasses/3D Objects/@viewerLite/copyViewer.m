function copyViewer(v,vref,name)
         if nargin < 3, name = 'result';end
         syncCamera(vref,v);
         v.BackgroundColor = vref.BackgroundColor;
         v.SceneLightVisible = vref.SceneLightVisible;
         v.SceneLightLinked = vref.SceneLightLinked;
         v.SceneLightPosition = vref.SceneLightPosition;
         set(v.Figure,'Name',name);
         set(v.Figure,'Position',get(vref.Figure,'Position'));
end