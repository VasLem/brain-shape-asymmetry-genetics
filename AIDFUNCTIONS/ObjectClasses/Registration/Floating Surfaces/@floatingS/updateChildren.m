function updateChildren(obj,action,varargin)
         switch obj.SurfaceClass
             case {'meshObj' 'LMObj'}
                 updateChildren(obj.Surface,action,varargin{:});
             otherwise
                 return;
         end
end