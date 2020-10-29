function centerVertices(obj)
% setting the vertices around their gravity point
  gravity = mean(obj.Vertices,2);
  obj.Vertices = obj.Vertices - repmat(gravity,1,obj.nrV);
  if ~isempty(obj.PoseLM), obj.PoseLM.Vertices = obj.PoseLM.Vertices-repmat(gravity,1,obj.PoseLM.nrV); end
      
end