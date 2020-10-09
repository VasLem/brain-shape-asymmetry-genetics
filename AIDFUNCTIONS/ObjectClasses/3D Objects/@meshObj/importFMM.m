function importFMM(obj,filename,varargin)
         [path, type] = readVarargin(varargin{:});
         cd(path);
         load(filename);
         switch type
             case 'FMM Original'
                 mesh = scan.mesh;
                 mesh.Mapped = false;
             case 'FMM Mapped';
                 mesh = scan.reg_mesh;
                 mesh.Mapped = true;
         end
         %clear scan;
         obj.Vertices = mesh.Location;
         obj.Faces = mesh.Tri;
         obj.ColorMode = 'Single';
         if isfield(mesh,'Color')
              if isfield(mesh.Color,'Value')
                 obj.TextureColor = mesh.Color.Value;
                 try
                     obj.UV = mesh.Color.uv;
                     obj.TextureMap = mesh.Color.im;
                 catch
                 end
                 obj.ColorMode = 'Texture';
              end
         end
         if isfield(mesh,'rbf')
            obj.RBF = fastRBF;
            obj.RBF.P = mesh.rbf;
         end
         LM = LMObj;
         LM.Vertices = scan.pose_LM;
         obj.PoseLM = LM;
end

function [path, type] = readVarargin(varargin)
  Input = find(strcmp(varargin,'Path'));
  if isempty(Input)
      path = pwd;
  else
      path = varargin{Input+1};
  end
  Input = find(strcmp(varargin,'Type'));
  if isempty(Input)
      type = 'FMM Original';
  else
      type = varargin{Input+1};
  end
end