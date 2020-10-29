function setCurrentMeshHighlight(obj)
  switch obj.Status
      case 'Ready'
         switch obj.SelectionMode
             case {'landmarks' 'areas' 'areasconnected' 'areasradius'}
                 if ~isempty(obj.CurrentMesh), obj.CurrentMesh.Border.SingleColor = [0.5 0 0]; end
             otherwise
                 meshes = obj.MeshChildren;
                 if isempty(meshes), return; end
                 for m=1:1:length(meshes)
                     meshes{m}.Border.SingleColor = [0.5 0.5 0.5];
                 end
         end
      case 'Busy'
      otherwise
  end
end