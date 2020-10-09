classdef borderObj < patchObj
    properties
        % Border properties
        VerticesIndex = [];
        RBF = [];
        ViewMode = 'Points';
        IndexedColor = [];
    end
    properties (Dependent = true)
        Vertices;
        Faces;
        TextureColor;
        %IndexedColor;
        Distance2Border;
        EdgeDistance2Border;
    end     
    methods %Constructor
        function obj = borderObj(varargin)
          obj = obj@patchObj(varargin{:});   % call patchObj constructor  
          if nargin > 0
           Input = find(strcmp(varargin, 'VerticesIndex'));if~isempty(Input), obj.VerticesIndex = varargin{Input+1}; end
           Input = find(strcmp(varargin, 'RBF'));if~isempty(Input), obj.RBF = varargin{Input+1}; end
           Input = find(strcmp(varargin, 'Distance2Border'));if~isempty(Input), obj.Distance2Border = varargin{Input+1}; end
           Input = find(strcmp(varargin, 'EdgeDistance2Border'));if~isempty(Input), obj.EdgeDistance2Border = varargin{Input+1}; end
          end
          obj.MarkerSize = 5;
          if validParent(obj) && strcmp(class(obj.Parent),'meshObj'), 
             obj.Axes = obj.Parent.Axes; 
             switch obj.Parent.Selected
              case true
                   obj.SingleColor = [0.5 0 0];
              case false
                   obj.SingleColor = [0.5 0.5 0.5];
             end
          end            
        end % LMObj Constructor
    end  
    methods % Special Setting and Getting
        function faces = get.Faces(obj)
            if isempty(obj.Vertices), faces = []; return; end
            faces = [(1:1:size(obj.Vertices,2));nan*ones(1,size(obj.Vertices,2))];
        end
        function obj = set.Faces(obj,faces) %#ok<INUSD>
            % do nothing
        end
        function vertices = get.Vertices(obj)
            if ~validParent(obj), vertices = []; return; end
            vertices = obj.Parent.Vertices(:,obj.VerticesIndex);
        end
        function obj = set.Vertices(obj,Loc)
            if ~(size(Loc,1)==3), Loc = Loc'; end
            if isempty(Loc), return; end% can't be done
            obj.Parent.Vertices(:,obj.VerticesIndex) = Loc;
            checkVerticesColorUpdate(obj,'update vertices');
        end
        function obj = set.VerticesIndex(obj,index)
                 if isempty(index) % special care
                    deletePatch(obj); 
                    obj.VerticesIndex = [];
                    %obj.TextureColor = [];
                    obj.IndexedColor = [];
                    return;
                 end
                 obj.RBF = [];% if index changes, then the current rbf is not valid anymore
                 obj.VerticesIndex = index;
                 checkVerticesColorUpdate(obj,'update vertices');
                 if ~validParent(obj), return; end
                 %if ~isempty(obj.Parent.TextureColor),obj.TextureColor = obj.Parent.TextureColor(:,index);end
                 %if ~isempty(obj.Parent.IndexedColor),obj.IndexedColor = obj.Parent.IndexedColor(:,index);end               
        end      
        function obj = set.ViewMode(obj,mode)
                %warning('Cannot change ViewMode for Border objects: Action Ignored');
                return; %dummy
        end
        function ViewMode = get.ViewMode(obj)
            ViewMode = obj.ViewMode;
        end
        function distances = get.Distance2Border(obj)
            distances = [];
            if isempty(obj.VerticesIndex)||~validParent(obj), return; end
            distances = intraDistances(obj.Parent,'VertexIndex',obj.VerticesIndex,'Type','triangle');
        end
        function obj = set.Distance2Border(obj,dist)
            % do nothing
        end
        function distances = get.EdgeDistance2Border(obj)
            distances = [];
            if isempty(obj.VerticesIndex)||~validParent(obj), return; end
            distances = intraDistances(obj.Parent,'VertexIndex',obj.VerticesIndex,'Type','edge');
        end
        function obj = set.EdgeDistance2Border(obj,dist)
            % do nothing
        end
        function obj = set.RBF(obj,rbf)
                 if ~isempty(obj.RBF); delete(obj.RBF); end
                 obj.RBF = rbf; 
        end
        function rbf = get.RBF(obj)
            rbf = obj.RBF;
            if ~mySuperClass.isH(rbf), rbf = []; return; end
        end
        function out = get.TextureColor(obj)
            if ~validParent(obj), out = []; return; end
            if isempty(obj.Parent.TextureColor), out = []; return; end
            out = obj.Parent.TextureColor(:,obj.VerticesIndex);
        end
        function obj = set.TextureColor(obj,in)
            % Do nothing;
        end
    end
    methods % Delete
%         function delete(obj)          
%            deletePatch(obj);
%         end
    end
end % classdef

