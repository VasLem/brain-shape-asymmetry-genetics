function zoomCamera(obj,xy,q)
         if isempty(q), q = max(-.9, min(.9, sum(xy)/70)); q = 1+q; end
        % hueristic avoids small view angles which will crash on solaris
         MIN_VIEW_ANGLE = .001;
         MAX_VIEW_ANGLE = 75;
         vaOld = camva(obj.RenderAxes);
         camzoom(obj.RenderAxes,q);
	     va = camva(obj.RenderAxes);
        %If the act of zooming puts us at an extreme, back the zoom out
         if ~((q>1 || va<MAX_VIEW_ANGLE) && (va>MIN_VIEW_ANGLE))
            set(obj.RenderAxes,'CameraViewAngle',vaOld);
         end
%          if obj.SceneLightLinked
%              set(obj.SceneLight,'Position',get(obj.RenderAxes,'CameraPosition'));
%          end
         drawnow expose;
end
