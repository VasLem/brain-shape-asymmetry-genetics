function panCamera(obj,xy)
         xy = -xy;
         xy = xy*camva(obj.RenderAxes)/500;
         campan(obj.RenderAxes,xy(1),xy(2),'none')%unconstrained panning
         if obj.SceneLightLinked
             set(obj.SceneLight,'Position',get(obj.RenderAxes,'CameraPosition'));
         end
         drawnow expose;
end