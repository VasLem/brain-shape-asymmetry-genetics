function rotateCamera(obj,xy)
         xy = -xy;
         camorbit(obj.RenderAxes,xy(1),xy(2),'none')%unconstrained rotation
         if obj.SceneLightLinked
             set(obj.SceneLight,'Position',get(obj.RenderAxes,'CameraPosition'));
         end
         drawnow expose;
end