classdef LMObj < patchObj
    properties
        % LM properties
        Vertices = [];
        TextureColor = [];
        IndexedColor = [];
        ViewMode = 'Points';
        Names = {};
        Vi = [];% landmarks on a mesh with 3 closest vertex (indices)
        Fi = [];% landmarks on a mesh with the index of the enclosing face
        Bary = [];
    end
    properties (Dependent = true)
        Faces;% LM objects don't have explicit Faces
    end
    methods %Constructor
        function obj = LMObj(varargin)
          obj = obj@patchObj(varargin{:});   % call patchObj constructor
          if nargin>0
           Input = find(strcmp(varargin, 'Names'));if~isempty(Input), obj.Names = varargin{Input+1}; end
           Input = find(strcmp(varargin, 'Vi'));if~isempty(Input), obj.Vi = varargin{Input+1}; end
           Input = find(strcmp(varargin, 'Fi'));if~isempty(Input), obj.Fi = varargin{Input+1}; end
          end          
          obj.MarkerSize = 20;
          obj.SingleColor = [0 0 1];
          obj.Type = 'LMObj';
        end % LMObj Constructor
    end  
    methods % Special Setting and Getting
        function obj = set.Faces(obj,faces) %#ok<INUSD>
                 % Do nothing
        end
        function faces = get.Faces(obj)
            if isempty(obj.Vertices), faces = []; return; end
            faces = [(1:1:size(obj.Vertices,2));nan*ones(1,size(obj.Vertices,2))];
        end
        function obj = set.Vertices(obj,Loc)
            if ~(size(Loc,1)==3), Loc = Loc'; end
            obj.Vertices = Loc;
            if isempty(Loc)% special care
                deletePatch(obj);
            end
            checkVerticesColorUpdate(obj,'update vertices');
        end
        function obj = set.ViewMode(obj,mode) %#ok<INUSD>
                % warning('Cannot change ViewMode for LM objects: action ignored');
                % warning('Cannot change ViewMode for LM objects');
        end
        function obj = set.TextureColor(obj,Color)
                 if ~(size(Color,1)==3), Color = Color'; end
                 if ~(size(Color,2)==size(obj.Vertices,2)), return; end
                 obj.TextureColor = Color;
                 checkVerticesColorUpdate(obj,'update texture');
        end
        function obj = set.IndexedColor(obj,Color)
                 if ~(size(Color,1)==1), Color = Color'; end
                 obj.IndexedColor = Color;
                 checkVerticesColorUpdate(obj,'update index');
        end
    end
    methods % barycentric coordinate functions
        function obj = cart2bary(obj,mesh)
                 obj.Bary = cart2bary(mesh.Vertices,mesh.Faces,obj.Vertices,obj.Fi); 
        end
        function obj = bary2cart(obj,mesh)
                 obj.Vertices = bary2cart(mesh.Vertices,mesh.Faces,obj.Bary,obj.Fi);
        end
        function obj = baryValue(obj,mesh)
                 if isempty(obj.Bary), cart2bary(obj,mesh); end
                 obj.Value = zeros(1,length(obj.Fi));
                 for i=1:1:length(obj.Fi)
                    f = obj.Fi(i);
                    fp = mesh.Value(:,mesh.Faces(:,f));
                    obj.Value(i) = obj.Bary(1,i)*fp(:,1)+obj.Bary(2,i)*fp(:,2)+obj.Bary(3,i)*fp(:,3);
                 end
        end
    end
    methods (Static = true)
        function obj = import(filename, varargin)
                 
        end
    end
end % classdef

