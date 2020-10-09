classdef patchObj < superClass
    % This is a class to create 3D Objects, which are linked with a graphical
    % patch object. The idea is that property changes of the 3D object reflect 
    % visually through the underlying patch obj. Furthermore, the underlying
    % patch can be deleted or destroyed when closing the figure containing
    % the patch object and further processing is still possible.
    % As is the patch handle stored in the object when visualized, so is
    % the object handle stored in the patch UserData field. This enables
    % visual interaction and modification of the object.
    %
    % It is a subclass of the handles hgsetget class, such that the objects
    % are handles themselves and such that the set get interface can be
    % used.
    %
    % This is the superclass for subclasses surfaceObj and lmObj, 
    %properties (Abstract=true)% are defined in subclasses
%     properties (GetAccess = private, SetAccess = private) 
%         Vertices;% Vector containing Vertices points or Vertices
%         Faces;% Vector containing the Vertices connections
%         ViewMode;% Defines how the patch is visualized (Solid, wireframe, points, ...)\
%         TextureColor;% Defines the texture color on the surface
%         IndexedColor = [];% Indexed Color value per Vertices point
%     end
    properties (Abstract = true) 
        Vertices;% Vector containing Vertices points or Vertices
        Faces;% Vector containing the Vertices connections
        ViewMode;% Defines how the patch is visualized (Solid, wireframe, points, ...)\
        TextureColor;% Defines the texture color on the surface
        IndexedColor;% Indexed Color value per Vertices point
    end
    properties
        % 3D object mesh and color properties
            %TextureColor = [];% Different RGB values per Vertices point
            %IndexedColor = [];% Indexed Color value per Vertices point
            SingleColor = [0 0 0];% Single RGB color per Vertices point
            Distance = [];% Extra Field for Distance information
            Indices = [];% Field to store active Vertices of 3D object, commonly used in registration objects
        % Patch & Rendering properties
            ColorMode = 'Single';% Defines the coloring of the patch
            Alpha = 1;% Defines the transparency of the patch
            Material = 'Default';% Defines how the patch visually behaves when Light is in the scene
            LightMode = 'flat';% Defines the reflection quality
            MarkerSize = 5;
            %ph = [];% This is the underlying patch handle     
            UpdatePatch = true;
        % General properties
            Tag = ''; % Tag can be used as a naming property
            Parent = [];% Parent of the obj
            UserData = [];% Userdefined data storage
            Selected = false;% Selected, to determine whether object is active or not e.g.
            %Visible =  false;% Show whether the object must be visible or not
    end
    properties (Dependent = true)
        nrV; %number of vertices;
        nrF; %number of faces;
        Gradient; % is taken care for in the patch object
        Location; % is linked with Vertices (fastRBF interface)
        Tri% Linked with faces (fastRBF interface)
        Value;% Linked with indexedColor
        RGB;% linked with TextureColor      
    end
    properties (Transient = true)
            %Selected = false;
            Visible = false;
            ph = [];% This is the underlying patch handle
            Axes = [];% A patch handle requires an axes object for creation
    end
    methods %Constructor
        function obj = patchObj(varargin)
          obj = obj@superClass(varargin{:});
          obj.SingleColor = 0.2 + (1-0.2).*rand(1,3);
        end
    end  
    methods % Special Setting and Getting
        function obj = set.Visible(obj,show)
            if~islogical(show), error('Visible must be true or false');end
            obj.Visible = show;                      
            if show
               if isempty(obj.ph)||~ishandle(obj.ph), createPatch(obj); end
            end
            setPatch(obj,'Visible');
            updateChildren(obj,'Visible Change');
        end
        function out = get.Visible(obj)
            if isempty(obj.ph)||~ishandle(obj.ph)||~obj.UpdatePatch, out = false; return; end
            out = obj.Visible;
        end
        function grad = get.Gradient(obj)
                 if strcmp(obj.Type,'image3D'), grad = [];return; end
                 if isempty(obj.Vertices), grad = []; return; end
                 if isempty(obj.Faces), grad = []; return; end
                 if ~isempty(find(isnan(obj.Faces), 1)), grad = []; return; end
                 tmp.vertices = obj.Vertices';tmp.faces = obj.Faces';
                 grad = patchnormals(tmp)';
                 grad = -1.*grad./repmat(sqrt(sum(abs(grad).^2,1)),size(grad,1),1);
        end
        function obj = set.SingleColor(obj,Color)
                 if ~(numel(Color)==3), error('Invalid single Color argument'); end
                 if ~(size(Color,2)==3), Color = Color'; end
                 obj.SingleColor = Color;
                 setPatch(obj,'ColorMode');
        end
        function obj = set.ColorMode(obj,mode)
                 switch mode
                     case 'Texture'
                         if ~(size(obj.TextureColor,2)==size(obj.Vertices,2)), return; end
                     case {'Indexed' 'Values'}
                         if ~(length(obj.IndexedColor)==size(obj.Vertices,2)), return; end
                     case 'Single'
                     otherwise
                      error('ColorMode must be Texture, Indexed (Values) or Single');   
                 end
                 obj.ColorMode = mode;
                 setPatch(obj,'ColorMode');
        end
        function obj = set.Alpha(obj,val)
                 if (val<0) || (val>1)
                    error('Alpha Value must be in the range of [0 1]');
                 end
                 obj.Alpha = val;
                 setPatch(obj,'Transparency');
        end
        function obj = set.Material(obj,mode)
                 if ~(strcmpi(mode,'Facial') ||... 
                      strcmpi(mode,'Dull') ||... 
                      strcmpi(mode,'Shiny') ||... 
                      strcmpi(mode,'Metal') ||...
                      strcmpi(mode,'Image3D') ||...
                      strcmpi(mode,'Default'))
                      error('Material must be Facial, Dull, Shiny, Metal or Default');
                 end
                 obj.Material = mode;
                 setPatch(obj,'Material');
        end
        function obj = set.LightMode(obj,mode)
                 if ~(strcmpi(mode,'none') ||... 
                      strcmpi(mode,'flat') ||... 
                      strcmpi(mode,'gouraud') ||... 
                      strcmpi(mode,'phong'))
                      error('Lighting must be none, flat, gouraud or phong');
                 end
                 obj.LightMode = mode;
                 setPatch(obj,'LightMode');
        end
        function obj = set.Axes(obj,ax)
            if isempty(ax), obj.Axes = []; return; end
            if ~ishandle(ax), obj.Axes = []; return; end
            if ~strcmp(get(ax,'Type'),'axes'),obj.Axes = []; return; end
            obj.Axes = ax;
            setPatch(obj,'Axes');
            updateChildren(obj,'Axe Change');
        end
        function obj = set.MarkerSize(obj,ms)
            obj.MarkerSize = ms;
            setPatch(obj,'MarkerSize');
        end
        function obj = set.UpdatePatch(obj,up)
            switch up
                case true
                    obj.UpdatePatch = true;
                case false
                    obj.UpdatePatch = false;
                    return;
                otherwise
                    return;
            end
            setPatch(obj,'All');
        end
        function obj = set.Selected(obj,log)
                 if ~islogical(log), return; end
                 obj.Selected = log;
                 updateChildren(obj,'Selected Change');
        end
        function nrV = get.nrV(obj)
                 nrV = size(obj.Vertices,2);
        end
        function nrF = get.nrF(obj)
                 nrF = size(obj.Faces,2);
        end
        function loc = get.Location(obj)
            loc = obj.Vertices;
        end
        function obj = set.Location(obj,loc)
            obj.Vertices = loc;
        end
        function tri = get.Tri(obj)
            tri = obj.Faces;
        end
        function obj = set.Tri(obj,tri)
            obj.Faces = tri;
        end
        function val = get.Value(obj)
            val = obj.IndexedColor;
        end
        function obj = set.Value(obj,val)
            obj.IndexedColor = val;
        end
        function rgb = get.RGB(obj)
            rgb = obj.TextureColor;
        end
        function obj = set.RGB(obj)
            obj.TextureColor = rgb;
        end
        function obj = set.Distance(obj,dist)
                 if numel(dist) == 1
                    dist = repmat(dist,1,obj.nrV);
                 end
                 obj.Distance = dist;
        end
        function dist = get.Distance(obj)
                 dist = obj.Distance;
                 %obj.nrV
                 if isempty(dist)
                    dist = zeros(1,obj.nrV);
                 end
        end
        function ind = get.Indices(obj)
            ind = obj.Indices;
            if isempty(ind), ind = (1:obj.nrV); end
        end
    end
    methods % Interface functions
        function delete(obj)
           deletePatch(obj);
           delete@superClass(obj); 
        end
        function cobj = fastClone(obj)
            % only clones surface information
            % Adviced to use when applying (implementing) functions which change ONLY
            % limited information (ex. TM (Transformation) class)
            cobj = eval(class(obj));
            fastCopy(obj,cobj);
        end
        function fastCopy(obj,cobj)
            % only copying surface information
            % Adviced to use when applying (implementing) functions which change ONLY
            % limited information (ex. TM (Transformation) class)
            if ~strcmp(obj.Type,cobj.Type), error('Copy: objects must be of same type'); end
            cobj.Vertices = obj.Vertices;
            cobj.Faces = obj.Faces;
            cobj.TextureColor = obj.TextureColor;
            cobj.IndexedColor = obj.IndexedColor;
            cobj.Distance = obj.Distance;
        end        
        function dist = vDistances(obj,obj2)
            % Distances between corresponding vertices
            if ~(obj.nrV==obj2.nrV), error('Number of vertices for both objects need to be the same to take difference'); end
            dist = sqrt(sum((obj.Vertices-obj2.Vertices).^2));
        end
        function diff = vDifferences(obj,obj2)
            % Distances between corresponding vertices
            if ~(obj.nrV==obj2.nrV), error('Number of vertices for both objects need to be the same to take difference'); end
            diff = obj.Vertices-obj2.Vertices;
        end
        function rmedse = vRmedse(obj,obj2)
            % RMedSe between corresponding vertices
                 rmedse = sqrt(med(vDistances(obj,obj2).^2));
        end
        function rmse = vRmse(obj,obj2)
            % RMedSe between corresponding vertices
                 rmse = sqrt(mean(vDistances(obj,obj2).^2));
        end
        
    end
    methods % validity children
        function out = validChild(obj,child) %#ok<INUSL>
            if isempty(child), out = false; return; end
            if ~isa(child,'handle'), out = false; return; end
            if ~isvalid(child), out = false; return; end
            if ~(child.Parent == obj), out = false; return; end
            out = true;
        end
        function out = validParent(obj) %#ok<INUSL>
            if isempty(obj.Parent), out = false; return; end
            if ~isa(obj.Parent,'handle'), out = false; return; end
            if ~isvalid(obj.Parent), out = false; return; end
            out = true;
        end 
    end
    methods (Static = true)
        function out = sum(objlist,field,w)
            %weighted sum of object Fields
            nrO = length(objlist);
            if nrO == 1, out = fastClone(objlist{1}); return; end
            out = fastClone(objlist{1});
            if nargin < 4,w = ones(nrO,out.nrV);end
            if size(w,2)==1, w = repmat(w,1,out.nrV);end
            sumW = sum(w,1);
            switch field
                case 'Vertices'
                    out.Vertices = repmat(w(1,:),3,1).*out.Vertices;
                    for i=2:1:nrO
                        out.Vertices = out.Vertices + repmat(w(i,:),3,1).*objlist{i}.Vertices;
                    end
                    out.Vertices = out.Vertices./repmat(sumW,3,1);
                case 'TextureColor'
                    out.TextureColor =repmat(w(1,:),3,1).*out.TextureColor;
                    for i=2:1:nrO
                        out.TextureColor = out.TextureColor + repmat(w(i,:),3,1).*objlist{i}.TextureColor;
                    end
                    out.TextureColor = out.TextureColor./repmat(sumW,3,1);
                case 'Distance'
                    out.Distance = w(1,:).*out.Distance;
                    for i=2:1:nrO
                        out.Distance = out.Distance + w(i,:).*objlist{i}.Distance;
                    end
                    out.Distance = out.Distance./sumW;
                case 'Value'
                    out.Value = w(1,:).*out.Value;
                    for i=2:1:nrO
                        out.Value = out.Value + w(i,:).*objlist{i}.Value;
                    end
                    out.Value = out.Value./sumW;
                otherwise
                    error('Invalid Field argument for summing two 3D object')
            end
        end
    end
end % classdef

%         function obj = set.IndexedColor(obj,Color)
%                  if ~(size(Color,1)==1), Color = Color'; end
%                  obj.IndexedColor = Color;
%                  checkVerticesColorUpdate(obj,'update index');
%         end

%         function grad = get.Gradient(obj)
%                  %grad = [];
%                  disp('retrieving gradient');
%                  if strcmp(obj.Type,'image3D'), grad = [];return; end
%                  if isempty(obj.Vertices), grad = []; return; end
%                  if isempty(obj.Faces), grad = []; return; end
%                  if ~isempty(find(isnan(obj.Faces), 1)), grad = []; return; end
%                  if isempty(obj.ph)||~ishandle(obj.ph)
%                     % in the background
%                     oldAxes = obj.Axes;
%                     f = figure('Visible','off');
%                     obj.Axes = gca;
%                     createPatch(obj);
%                     obj.Axes = gca;
%                     l = camlight;
%                     grad = get(obj.ph,'VertexNormals')';
%                     delete(l);
%                     delete(obj.ph);
%                     delete(f);
%                     obj.Axes = oldAxes;
%                  else
%                     grad = get(obj.ph,'VertexNormals')';
%                     if isempty(grad)% turn on the light
%                        axes(obj.Axes);l = camlight;
%                        grad = get(obj.ph,'VertexNormals')';
%                        delete(l);
%                     end
%                  end
%                  grad = grad./repmat(sqrt(sum(abs(grad).^2,1)),size(grad,1),1);
%         end






%           if nargin > 0 
%             Input = find(strcmp(varargin, 'Vertices'));if ~isempty(Input), obj.Vertices = varargin{Input+1}; end
%             Input = find(strcmp(varargin, 'Location'));if ~isempty(Input), obj.Location = varargin{Input+1}; end 
%             Input = find(strcmp(varargin, 'Faces'));if ~isempty(Input), obj.Faces = varargin{Input+1}; end
%             Input = find(strcmp(varargin, 'Tri'));if ~isempty(Input), obj.Tri = varargin{Input+1}; end            
%             Input = find(strcmp(varargin, 'Tag'));if ~isempty(Input), obj.Tag = varargin{Input+1}; end
%             Input = find(strcmp(varargin, 'Parent'));if ~isempty(Input), obj.Parent = varargin{Input+1}; end
%             Input = find(strcmp(varargin, 'UserData'));if ~isempty(Input), obj.UserData = varargin{Input+1}; end
%             Input = find(strcmp(varargin, 'Selected'));if ~isempty(Input), obj.Selected = varargin{Input+1}; end
%             Input = find(strcmp(varargin, 'TextureColor'));if ~isempty(Input), obj.TextureColor = varargin{Input+1}; end
%             Input = find(strcmp(varargin, 'RGB'));if ~isempty(Input), obj.RGB = varargin{Input+1}; end
%             Input = find(strcmp(varargin, 'IndexedColor'));if ~isempty(Input), obj.IndexedColor = varargin{Input+1}; end
%             Input = find(strcmp(varargin, 'Value'));if ~isempty(Input), obj.Value = varargin{Input+1}; end
%             Input = find(strcmp(varargin, 'SingleColor'));if ~isempty(Input), obj.SingleColor = varargin{Input+1}; end
%             Input = find(strcmp(varargin, 'Distance'));if ~isempty(Input), obj.Distance = varargin{Input+1}; end
%             Input = find(strcmp(varargin, 'Indices'));if ~isempty(Input), obj.Indices = varargin{Input+1}; end
%             %Input = find(strcmp(varargin, 'ph'));if ~isempty(Input), obj.ph = varargin{Input+1}; end
%             Input = find(strcmp(varargin, 'ColorMode'));if ~isempty(Input), obj.ColorMode = varargin{Input+1}; end
%             Input = find(strcmp(varargin, 'ViewMode'));if ~isempty(Input), obj.ViewMode = varargin{Input+1}; end
%             Input = find(strcmp(varargin, 'Alpha'));if ~isempty(Input), obj.Alpha = varargin{Input+1}; end
%             Input = find(strcmp(varargin, 'Material'));if ~isempty(Input), obj.Material = varargin{Input+1}; end
%             Input = find(strcmp(varargin, 'LightMode'));if ~isempty(Input), obj.LightMode = varargin{Input+1}; end
%             Input = find(strcmp(varargin, 'Axes'));if ~isempty(Input), obj.Axes = varargin{Input+1}; end
%             Input = find(strcmp(varargin, 'Visible'));if ~isempty(Input), obj.Visible = varargin{Input+1}; end
%           end

%         function obj = set.TextureColor(obj,Color)
%                  if ~(size(Color,1)==3), Color = Color'; end
%                  if ~(size(Color,2)==size(obj.Vertices,2)), return; end
%                  obj.TextureColor = Color;
%                  checkVerticesColorUpdate(obj,'update texture');
%         end

%         function cobj = clone(obj)
%             cobj = clone@mySuperClass(obj);
%             cobj.Axes = obj.Axes;
%             %cobj.Visible = obj.Visible;           
%         end

%         function copy(obj,cobj)
%             copy@mySuperClass(obj,cobj);
%             cobj.Vertices = obj.Vertices;
%             cobj.Faces = obj.Faces;
%             cobj.ViewMode = obj.ViewMode;
%             cobj.TextureColor = obj.TextureColor;
%             cobj.IndexedColor = obj.IndexedColor;
%             cobj.Distance = obj.Distance;
%             cobj.ColorMode = obj.ColorMode;
%             cobj.Alpha = obj.Alpha;
%             cobj.Material = obj.Material;
%             cobj.LightMode = obj.LightMode;
%             cobj.MarkerSize = obj.MarkerSize;
%             cobj.UpdatePatch = obj.UpdatePatch;
%             cobj.Tag = obj.Tag;
%             cobj.Parent = obj.Parent;
%             cobj.UserData = obj.UserData;
%             cobj.Selected = obj.Selected;
%             cobj.Indices = obj.Indices;
%         end

%         function struc = obj2struc(obj)
%             struc.Type = obj.Type;
%             struc.Vertices = obj.Vertices;
%             %struc.Gradient = obj.Gradient;
%             struc.Faces = obj.Faces;
%             struc.ViewMode = obj.ViewMode;
%             struc.TextureColor = obj.TextureColor;
%             struc.IndexedColor = obj.IndexedColor;
%             struc.Distance = obj.Distance;
%             struc.Indices = obj.Indices;
%             struc.ColorMode = obj.ColorMode;
%             struc.Alpha = obj.Alpha;
%             struc.Material = obj.Material;
%             struc.LightMode = obj.LightMode;
%             struc.MarkerSize = obj.MarkerSize;
%             struc.UpdatePatch = obj.UpdatePatch;
%             struc.Tag = obj.Tag;
%             %struc.Parent = obj.Parent;
%             struc.UserData = obj.UserData;
%             struc.Selected = obj.Selected;
%             %struc.Visible = obj.Visible;
%         end
%         function obj = struc2obj(obj,struc)
%             obj.Type = struc.Type;
%             obj.Vertices = struc.Vertices;
%             obj.Faces = struc.Faces;
%             obj.ViewMode = struc.ViewMode;
%             obj.TextureColor = struc.TextureColor;
%             obj.IndexedColor = struc.IndexedColor;
%             obj.Distance = struc.Distance;
%             obj.ColorMode = struc.ColorMode;
%             obj.Alpha = struc.Alpha;
%             obj.Material = struc.Material;
%             obj.LightMode = struc.LightMode;
%             obj.MarkerSize = struc.MarkerSize;
%             obj.UpdatePatch = struc.UpdatePatch;
%             obj.Tag = struc.Tag;
%             %obj.Parent = struc.Parent;
%             obj.UserData = struc.UserData;
%             if isfield(struc,'Indices'), obj.Indices = struc.Indices; end% Remove later
%             obj.Selected = struc.Selected;
%             %obj.Visible = struc.Visible;           
%         end