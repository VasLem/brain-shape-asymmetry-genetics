function resetCamera(obj)
         camlookat(obj.RenderAxes);
%          if ~isempty(obj.CurrentMesh), camlookat(obj.CurrentMesh); end
         if obj.SceneLightLinked
             set(obj.SceneLight,'Position',get(obj.RenderAxes,'CameraPosition'));
         end
end