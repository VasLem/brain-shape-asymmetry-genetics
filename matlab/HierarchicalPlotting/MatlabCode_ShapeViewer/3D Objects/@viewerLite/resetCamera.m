function resetCamera(obj)
         camlookat(obj.RenderAxes);
         if obj.SceneLightLinked
             set(obj.SceneLight,'Position',get(obj.RenderAxes,'CameraPosition'));
         end
end