classdef viewerLite < hgsetget
    % this class is a class of 3D viewers
    properties
        % figure object
        Figure = [];
        % Displayed Objects
        DisplayObjects = [];
        % Mode status
        Mode = 'camera'% the current mode of the viewer (none, camera, light)
        SelectionMode = 'full'% the selection mode (none, landmarks, area)
        LandmarkSelection = [];
        AreaSelection = [];
        AreaSelectionRadius = [];
        SelectionSphere = [];
        % Light objects
        SceneLight = [];
        SceneLightRotMode = 'vertical';
        SceneLightLinked = false;
        % Axes object
        RenderAxes = [];
        % Toolbar objects
        Toolbar = [];
        % Context menu objects
        ContextMenu = [];
        % Filters
        FilterFigure = [];
        SmoothMode = 'combinatorial';
        SmoothRuns = 1;
        SubdivideMode = 'runs';
        SubdivideVal = 1;
        %SubdivideSize = 1;
        ReducePercentage = 5;
        % General
        MotionData = [];
        Action = 'none';
        ActiveKey = 'none';
        Parent = [];
        Selected = true;
        Status = 'Ready';
        Tag = '3D Viewer'; 
        NoCurrentMeshMsg = 'No active scan! To activate a scan Select AND Display it. When Multiple scans are selected and displayed the first is made active, it is adviced to select and display only one scan whitin this viewer'
        CalledProcess = 'none'
        Record = false;
        MovieFile = [];
        FrameCounter = 0;
        UserData= [];
        Application = [];
        Title = [];
        PrintValue = true;
    end
    properties (Dependent = true)
        BackgroundColor;
        Renderer;
        AxesVisible;
        AxesGrid;
        AxesBox;
        AxesWallColor;
        AxesXColor;
        AxesYColor;
        AxesZColor;
        CameraPosition;
        CameraTarget;
        CameraUpVector;
        CameraViewAngle;
        DataAspectRatio;
        SceneLightMode;
        SceneLightPosition;
        SceneLightColor;
        SceneLightVisible;
        Visible;
        MeshChildren;
        CurrentMesh;
        VoxelImageChildren;
        CurrentVoxelImage;
        CamProjection;
    end
    methods %Constructor
        function obj = viewerLite(varargin)
          if nargin>0
             Input = find(strcmp(varargin, 'Tag'));if ~isempty(Input), obj.Tag = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'Parent'));if ~isempty(Input), obj.Parent = varargin{Input+1}; end
          end
          constructViewer(obj);
          %set(obj.Toolbar.cam_mode,'State','on');
          %set(obj.Toolbar.full_mode,'State','on');
          %obj.SelectionSphere.Radius = 120;
          %obj.SelectionSphere.Center = [0; 0; 0];
          obj.BackgroundColor = [1 1 1];
          obj.SceneLightLinked = true;
          obj.SceneLightVisible = true;
          setMousePointer(obj);
          updateContextMenu(obj,'All');
          Input = find(strcmp(varargin, 'CalledProcess'));if ~isempty(Input), obj.CalledProcess = varargin{Input+1}; end 
        end % 3D viewer Constructor
    end  
    methods % Special Setting and Getting
        function color = get.BackgroundColor(obj)
            color = get(obj.Figure,'Color');
        end
        function obj = set.BackgroundColor(obj,color)
            if ~numel(color)==3, return; end
            if size(color,1)==3, color = color'; end
            set(obj.Figure,'Color',color);
        end
        function renderer = get.Renderer(obj)
            renderer = get(obj.Figure,'renderer');
        end
        function obj = set.Renderer(obj,renderer)
            if ~(strcmpi(renderer,'opengl') ||strcmpi(renderer,'zbuffer'))
                error('renderer must be opengl or zbuffer');
            end
            set(obj.Figure,'renderer',renderer);
            updateContextMenu(obj,'Renderer');
        end
        function axesvisible = get.AxesVisible(obj)
            switch get(obj.RenderAxes,'Visible')
                case 'on'
                    axesvisible = true;
                case 'off'
                    axesvisible = false;
            end
        end
        function obj = set.AxesVisible(obj,val)
            switch val
                case true
                    set(obj.RenderAxes,'Visible','on');
                case false
                    set(obj.RenderAxes,'Visible','off');
                otherwise
                    return;
            end
            updateContextMenu(obj,'Axes Visible')
        end
        function axesgrid = get.AxesGrid(obj)
            switch get(obj.RenderAxes,'XGrid')
                case 'on'
                    axesgrid = true;
                case 'off'
                    axesgrid = false;
            end
        end
        function obj = set.AxesGrid(obj,val)
            switch val
                case true
                    grid(obj.RenderAxes,'on');
                case false
                    grid(obj.RenderAxes,'off');
                otherwise
                    return;
            end
            updateContextMenu(obj,'Axes Grid');
        end
        function axesbox = get.AxesBox(obj)
            switch get(obj.RenderAxes,'Box')
                case 'on'
                    axesbox = true;
                case 'off'
                    axesbox = false;
            end
        end
        function obj = set.AxesBox(obj,val)
            switch val
                case true
                    set(obj.RenderAxes,'Box','on');
                case false
                    set(obj.RenderAxes,'Box','off');
                otherwise
                    return;
            end
            updateContextMenu(obj,'Axes Box');
        end
        function color = get.AxesWallColor(obj)
            color = get(obj.RenderAxes,'Color');
        end
        function obj = set.AxesWallColor(obj,color)
            if ~numel(color)==3, return; end
            if size(color,1)==3, color = color'; end
            set(obj.RenderAxes,'Color',color);
        end
        function color = get.AxesXColor(obj)
            color = get(obj.RenderAxes,'XColor');
        end
        function obj = set.AxesXColor(obj,color)
            if ~numel(color)==3, return; end
            if size(color,1)==3, color = color'; end
            set(obj.RenderAxes,'XColor',color);
        end
        function color = get.AxesYColor(obj)
            color = get(obj.RenderAxes,'YColor');
        end
        function obj = set.AxesYColor(obj,color)
            if ~numel(color)==3, return; end
            if size(color,1)==3, color = color'; end
            set(obj.RenderAxes,'YColor',color);
        end
        function color = get.AxesZColor(obj)
            color = get(obj.RenderAxes,'ZColor');
        end
        function obj = set.CamProjection(obj,proj)
            if ~(strcmp(proj,'orthographic')||strcmp(proj,'perspective'))
               error('Projection be must either orthographic or perspective');
            end
            set(obj.RenderAxes,'Projection',proj);
            updateContextMenu(obj,'Camera Projection');
        end
        function proj = get.CamProjection(obj)            
            proj = get(obj.RenderAxes,'Projection');
        end
        function obj = set.AxesZColor(obj,color)
            if ~numel(color)==3, return; end
            if size(color,1)==3, color = color'; end
            set(obj.RenderAxes,'ZColor',color);
        end
        function show = get.Visible(obj)
            switch get(obj.Figure,'Visible')
                case 'on'
                    show = true; 
                case 'off'
                    show = false;
            end
        end
        function obj = set.Visible(obj,val)
            switch val
                case true
                    set(obj.Figure,'Visible','on');
                case false
                    set(obj.Figure,'Visible','off');
                otherwise
            end
        end
        function obj = set.Mode(obj,mode)
            if ~(strcmp(mode,'none')||strcmp(mode,'camera')||strcmp(mode,'light'))
               error('Viewer mode must be none, camera or light');
            end
            oldmode = obj.Mode;
            obj.Mode = mode;
            updateContextMenu(obj,'Mode');
            updateToolbar(obj,'Mode',oldmode);        
            setMousePointer(obj);
        end
        function obj = set.SelectionMode(obj,mode)
            if ~(strcmp(mode,'none')||strcmp(mode,'landmark')||strcmp(mode,'area')||strcmp(mode,'brush')||strcmp(mode,'fill')||strcmp(mode,'full'))
               error('Viewer selection mode must be none, landmarks or areas');
            end
            oldmode = obj.SelectionMode;
            obj.SelectionMode = mode;
            updateToolbar(obj,'SelectionMode',oldmode);
            updateContextMenu(obj,'SelectionMode');           
            setMousePointer(obj);
        end
        function obj = set.SmoothMode(obj,mode)
            if ~(strcmp(mode,'combinatorial')||strcmp(mode,'distance')||strcmp(mode,'conformal'))
               error('Viewer smooth mode must be combinatorial, distance or conformal');
            end
            obj.SmoothMode = mode;
            if ~ishandle(obj.FilterFigure), return; end
            set(findobj(obj.FilterFigure,'Tag',mode),'Value',1);
%             updateContextMenu(obj,'SmoothMode');
        end
        function obj = set.SmoothRuns(obj,runs)
                 if ~isinteger(runs)||(runs<=0)||(runs>10), return; end
                 obj.SmoothRuns = runs;
                 if ~ishandle(obj.FilterFigure), return; end
                 set(findobj(obj.FilterFigure,'Tag','SmoothRunsMenu'),'Value',runs);                 
        end
        function obj = set.SubdivideMode(obj,mode)
            if ~(strcmp(mode,'runs')||strcmp(mode,'size'))
               error('Viewer triangle subdivide mode must be runs or maximum size');
            end
            obj.SubdivideMode = mode;
            if ~ishandle(obj.FilterFigure), return; end
            set(findobj(obj.FilterFigure,'Tag',mode),'Value',1);
            switch mode
                case 'runs'
                 set(findobj(obj.FilterFigure,'Tag','SubdivideRunsMenu'),'Enable','on');
                 set(findobj(obj.FilterFigure,'Tag','SubdivideSizeMenu'),'Enable','off');
                case 'size'
                 set(findobj(obj.FilterFigure,'Tag','SubdivideRunsMenu'),'Enable','off');
                 set(findobj(obj.FilterFigure,'Tag','SubdivideSizeMenu'),'Enable','on');   
            end
        end
        function obj = set.SubdivideVal(obj,val)
                 if (val<=0)||(val>10), return; end
                 obj.SubdivideVal = val;
                 if ~ishandle(obj.FilterFigure), return; end
                 set(findobj(obj.FilterFigure,'Tag','SubdivideValMenu'),'Value',val);
        end
        function obj = set.ReducePercentage(obj,perc)
                 if ~~isinteger(perc)||(perc<5)||(perc>95)||~(mod(perc,5)==0), return; end
                 obj.ReducePercentage = perc;
                 if ~ishandle(obj.FilterFigure), return; end
                 set(findobj(obj.FilterFigure,'Tag','ReducePercentageMenu'),'Value',perc/5);
        end
        function mode = get.SceneLightMode(obj)
            mode = get(obj.SceneLight,'Style');
        end
        function obj = set.SceneLightMode(obj,mode)
            if ~(strcmp(mode,'infinite')||strcmp(mode,'local'))
               error('Viewer scene light mode must be infinite or local');
            end
            set(obj.SceneLight,'Style',mode);
        end
        function obj = set.SceneLightRotMode(obj,mode)
            if ~(strcmp(mode,'vertical')||strcmp(mode,'horizontal'))
               error('Viewer scene light rotation mode must be vertical or horizontal');
            end
            obj.SceneLightRotMode = mode;
        end
        function obj = set.Tag(obj,tag)
            if ~ischar(tag), return; end
            obj.Tag = tag;
            if ishandle(obj.Figure), set(obj.Figure,'Name',tag); end
            if ishandle(obj.FilterFigure), set(obj.FilterFigure,'Name', ['Filter Settings ' tag]); end
        end
        function pos = get.SceneLightPosition(obj)
            pos = get(obj.SceneLight,'Position');
        end
        function obj = set.SceneLightPosition(obj,pos)
            set(obj.SceneLight,'Position',pos);
        end
        function color = get.SceneLightColor(obj)
            color = get(obj.SceneLight,'Color');
        end
        function obj = set.SceneLightColor(obj,color)
            if ~numel(color)==3, return; end
            if size(color,1)==3, color = color'; end
            set(obj.SceneLight,'Color',color);
        end
        function scenelightvisible = get.SceneLightVisible(obj)
            switch get(obj.SceneLight,'Visible')
                case 'on'
                    scenelightvisible = true;
                case 'off'
                    scenelightvisible = false;
            end
        end
        function obj = set.SceneLightVisible(obj,val)
            switch val
                case true
                    set(obj.SceneLight,'Visible','on');
                case false
                    set(obj.SceneLight,'Visible','off');
                otherwise
                    return;
            end
            updateContextMenu(obj,'Scene Light');
        end
        function obj = set.SceneLightLinked(obj,log)
                 if ~islogical(log), return; end
                 obj.SceneLightLinked = log;
                 updateContextMenu(obj,'Scene Light Linked');
        end            
        function pos = get.CameraPosition(obj)
            pos = get(obj.RenderAxes,'cameraposition');
        end
        function obj = set.CameraPosition(obj,pos)
            set(obj.RenderAxes,'cameraposition',pos);
        end
        function target = get.CameraTarget(obj)
            target = get(obj.RenderAxes,'cameratarget');
        end
        function obj = set.CameraTarget(obj,target)
            set(obj.RenderAxes,'cameratarget',target);
        end
        function out = get.CameraViewAngle(obj)
            out = get(obj.RenderAxes,'CameraViewAngle');
        end
        function obj = set.CameraViewAngle(obj,in)
            set(obj.RenderAxes,'CameraViewAngle',in);
        end
        function daspr = get.DataAspectRatio(obj)
            daspr = get(obj.RenderAxes,'dataaspectratio');
        end
        function obj = set.DataAspectRatio(obj,daspr)
            set(obj.RenderAxes,'DataAspectRatio',daspr);
        end
        function upv = get.CameraUpVector(obj)
            upv = get(obj.RenderAxes,'cameraupvector');
        end
        function obj = set.CameraUpVector(obj,upv)
            set(obj.RenderAxes,'cameraupvector',upv);
        end
        function children = get.MeshChildren(obj)
            patches = findobj(get(obj.RenderAxes,'Children'),'flat','Type','patch');
            children = {};
            if isempty(patches), return; end
            counter = 0;
            for p=1:1:length(patches)
                UD = get(patches(p),'UserData');
                if strcmp(class(UD),'meshObj')%||strcmp(class(UD),'image3D')
                   counter = counter + 1;
                   children{counter} = UD; %#ok<AGROW>
                end
            end
        end
        function obj = set.MeshChildren(obj,children)
            if isempty(children), return; end
            for c=1:1:length(children)
                if strcmp(class(children{c}),'meshObj')||strcmp(class(children{c}),'image3D')
                   children{c}.Axes = obj.RenderAxes;
                end
            end
        end
        function currentmesh = get.CurrentMesh(obj)
                 currentmesh = [];
                 meshes = obj.MeshChildren;
                 if isempty(meshes), currentmesh = []; return; end
                 for m=1:1:length(meshes)
                     if (meshes{m}.Selected) && (meshes{m}.Visible)
                         currentmesh = meshes{m};
                         break;
                     end
                 end            
        end
        function obj = set.VoxelImageChildren(obj,children)
            if isempty(children), return; end
            for c=1:1:length(children)
                if strcmp(class(children{c}),'voxelImage')
                   children{c}.Axes = obj.RenderAxes;
                end
            end
        end
        function out = get.VoxelImageChildren(obj)
            patches = findobj(get(obj.RenderAxes,'Children'),'flat','Type','patch');
            out = {};
            if isempty(patches), return; end
            counter = 0;
            for p=1:1:length(patches)
                UD = get(patches(p),'UserData');
                if strcmp(class(UD),'voxelImage')
                   counter = counter + 1;
                   out{counter} = UD; %#ok<AGROW>
                end
            end
        end
        function out = get.CurrentVoxelImage(obj)
                 out = [];
                 images = obj.VoxelImageChildren;
                 if isempty(images), out = []; return; end
                 for m=1:1:length(images)
                     if (images{m}.Selected) && (images{m}.Visible)
                         out = images{m};
                         break;
                     end
                 end            
        end
        function obj = set.Status(obj,stat)
                 obj.Status = stat;
                 setMousePointer(obj);
        end
        function obj = set.Action(obj,action)
                 obj.Action = action;
                 setMousePointer(obj);
        end
        function obj = set.ActiveKey(obj,key)
                 obj.ActiveKey = key;
                 setMousePointer(obj);
        end
        function obj = set.CalledProcess(obj,process)
                 obj.CalledProcess = process;
                 switch process
                     case 'indicatePoseLM'
                         set(obj.Toolbar.landmark_mode,'State','on');
                     case 'Batch Processing'
                         %set(obj.Toolbar.brush_mode,'State','on');
                         set(obj.Toolbar.area_mode,'State','on');
                         set(obj.Toolbar.light_toggle,'State','on');
                         set(obj.Toolbar.link_toggle,'State','on');
                     otherwise
                         return;
                 end
                 set(obj.Figure,'Name',process);
                 
        end
        function out = get.Application(obj)
                 out = obj.Application;
                 if ~superClass.isH(out), out = []; end
        end
    end
    methods % validity children
        function out = validChild(obj,child) %#ok<INUSL>
            if isempty(child), out = false; return; end
            if ~isa(child,'handle'), out = false; return; end
            if ~isvalid(child), out = false; return; end
            out = true;
        end       
    end
    methods % Delete
        function delete(obj)
           if ishandle(obj.Figure), delete(obj.Figure);end
           if ishandle(obj.FilterFigure), delete(obj.FilterFigure);end
        end
        function deletePatchChildren(obj)
            patches = findobj(get(obj.RenderAxes,'Children'),'flat','Type','patch');
            if isempty(patches), return; end
            for p=1:1:length(patches)
                UD = get(patches(p),'UserData');
                try
                    delete(UD)
                catch
                end
            end
            
        end
        function hidePatchChildren(obj)
            patches = findobj(get(obj.RenderAxes,'Children'),'flat','Type','patch');
            if isempty(patches), return; end
            for p=1:1:length(patches)
                UD = get(patches(p),'UserData');
                try
                    UD.Visible = false;
                catch
                end
            end
            
        end
    end
end % classdef
