function [Features,Faces,Flags] = convert2MeshMonk(obj,Flags)
          Features = single([obj.Vertices',obj.Gradient']);
          Faces = uint32(obj.Faces'-1);
          if nargin<2, Flags = single(ones(obj.nrV,1));return;end
          Flags = single(Flags);   
end