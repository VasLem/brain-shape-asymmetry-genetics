classdef vectorField3D < patchObj
    properties
        % Vector properties
        StartPoints = [];
        EndPoints = [];
        ViewMode = 'Wireframe';
    end
    properties (Dependent = true)
        Faces;
        Vertices;
        TextureColor;
        IndexedColor;
        nrSP;
        nrEP;
        Direction;
        Length;
        Difference;
    end
    methods %Constructor
        function obj = vectorField3D(varargin)
          obj = obj@patchObj(varargin{:});   % call patchObj constructor
          obj.ColorMode = 'Indexed';
        end % LMObj Constructor
    end  
    methods % Special Setting and Getting
        function obj = set.Faces(obj,in) %#ok<INUSD>
                 % Do nothing
        end
        function obj = set.Vertices(obj,in)
            % Do nothing
        end
        function obj = set.TextureColor(obj,in)
                 % Do nothing
        end
        function obj = set.IndexedColor(obj,in)
        end 
        function obj = set.ViewMode(obj,in) %#ok<INUSD>
                % warning('Cannot change ViewMode for LM objects: action ignored');
                % warning('Cannot change ViewMode for LM objects');
        end
        function faces = get.Faces(obj)
            if isempty(obj.Vertices), faces = []; return; end
            faces = [(1:1:obj.nrEP);...
                     (obj.nrEP+1:1:2*obj.nrEP);...
                     nan*ones(1,obj.nrEP)];
        end    
        function out = get.Vertices(obj)
            if isempty(obj.EndPoints), out = []; return; end
            sp = obj.StartPoints;
            ep = obj.EndPoints;
            if ~size(sp,2)==obj.nrEP, warning('No equal number end and start points'), out = []; return; end
            out = [sp, ep];
        end
        function out = get.TextureColor(obj)
            if isempty(obj.Vertices), out = []; return; end
%             out = [[ones(1,obj.nrEP);zeros(1,obj.nrEP);zeros(1,obj.nrEP)] 0.5*ones(3,obj.nrEP)];
            out = [[zeros(1,obj.nrEP);zeros(1,obj.nrEP);ones(1,obj.nrEP)] [ones(1,obj.nrEP);zeros(1,obj.nrEP);zeros(1,obj.nrEP)]];
        end
        function out = get.IndexedColor(obj)
            if isempty(obj.Vertices), out = []; return; end
            %out = [zeros(1,obj.nrEP) obj.Length];
            out = [obj.Length obj.Length];
        end
        function out = get.nrSP(obj)
            out = size(obj.StartPoints,2);
        end
        function out = get.nrEP(obj)
            out = size(obj.EndPoints,2);
        end
        function out = get.Direction(obj)
            if isempty(obj.EndPoints), out = []; return; end
            out = obj.Difference./repmat(obj.Length,3,1);
        end
        function out = get.Length(obj)
            if isempty(obj.EndPoints), out = []; return; end
            out = sqrt(sum(obj.Difference.^2));
        end
        function out = get.Difference(obj)
            if isempty(obj.EndPoints), out = []; return; end
            out = obj.EndPoints - obj.StartPoints;            
        end
        function out = get.StartPoints(obj)
            if isempty(obj.EndPoints), out = []; return; end
            out = obj.StartPoints;
            if isempty(out), out = [0;0;0]; end
            if size(out,2) == 1, out = repmat(out,1,obj.nrEP); end
        end
        function obj = set.EndPoints(obj,in)
            obj.EndPoints = in;
            checkVerticesColorUpdate(obj,'update vertices');
        end
        function obj = set.StartPoints(obj,in)
            obj.StartPoints = in;
            checkVerticesColorUpdate(obj,'update vertices');
        end
    end
    methods % InterFace Functions
        function export(obj,filename,path)
            if nargin==3, cd(path);end
            mesh.Location = obj.Vertices;
            mesh.Tri = obj.Faces;
            fastrbf_export(mesh,filename,'obj');
        end
    end
end % classdef

