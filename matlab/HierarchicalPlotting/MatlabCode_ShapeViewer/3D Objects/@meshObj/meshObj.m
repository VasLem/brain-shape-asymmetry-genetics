classdef meshObj < patchObj
    properties
        % superclass abstract properties
            Vertices = [];
            Faces = [];
            TextureColor = [];
            IndexedColor = [];
            ViewMode = 'Solid';
        % mesh specific properties
            RBF = [];
            KDE = [];
            Border = [];
            TextureMap = [];
            UV = [];
            %UVClusters = [];
            Mapped = false;
            Reference = false;
        % LM children
            PoseLM = {};
            CustomLM = {};
        % Area children
            Area = {};
            LinkChildren = true;
    end
    properties (Dependent = true)
        Adjacency;
        nrCLM;
        nrA;
    end
    methods %Constructor
        function obj = meshObj(varargin)
          obj = obj@patchObj(varargin{:});   % call patchObj constructor
          if ~isempty(obj.TextureColor)
              obj.ColorMode = 'Texture';
          elseif ~isempty(obj.IndexedColor)
              obj.ColorMode = 'Indexed';
          end
          obj.LightMode = 'gouraud';
        end % LMObj Constructor
    end  
    methods % Special Setting and Getting
        function obj = set.Vertices(obj,Loc)
                 if ~(size(Loc,1)==3), Loc = Loc'; end
                 if size(Loc,2)==0, Loc = []; end
                 obj.Vertices = Loc;
                 checkVerticesColorUpdate(obj,'update vertices');
        end
        function obj = set.Faces(obj,Faces)
                 if ~(size(Faces,1)==3), Faces = Faces'; end
                 if size(Faces,2)==0||~(size(Faces,1)==3), Faces = []; end
                 obj.Faces = Faces;
                 setPatch(obj,'Faces');
        end
        function out = get.Faces(obj)
                 out = obj.Faces;
                 if ~isempty(out), return; end
                 if isempty(obj.Vertices), return; end
                 out = [(1:1:size(obj.Vertices,2));nan*ones(1,size(obj.Vertices,2))];
        end
        function obj = set.ViewMode(obj,mode)
                 if ~(strcmpi(mode,'Solid') ||... 
                      strcmpi(mode,'Wireframe') ||... 
                      strcmpi(mode,'Points') ||... 
                      strcmpi(mode,'Solid/Wireframe'))
                      error('ViewMode must be Solid, Wireframe, Solid/Wireframe or Points');
                 end
                 obj.ViewMode = mode;
                 setPatch(obj,'ViewMode');
        end
        function obj = set.RBF(obj,rbf)
                 if ~isempty(obj.RBF); delete(obj.RBF); end
                 obj.RBF = rbf; 
        end
        function rbf = get.RBF(obj)
            rbf = obj.RBF;
            if ~mySuperClass.isH(rbf), rbf = []; return; end
        end
        function obj = set.Border(obj,border)
                 if ~isempty(obj.Border), delete(obj.Border); end
                 if isempty(border), obj.Border = []; return; end
                 if ~strcmp(class(border),'borderObj'), error('no valid border object'); end                
                 obj.Border = border;
                 obj.Border.Parent = obj;
        end
        function out = get.Border(obj)
            out = obj.Border;
            if~mySuperClass.isH(out), out = []; return; end
            %if isempty(out.Vertices), out = []; return; end
        end
        function obj = set.TextureMap(obj,in)
                 if ~isempty(obj.TextureMap), delete(obj.TextureMap); end
                 if isempty(in), obj.TextureMap = []; return; end
                 switch class(in)
                     case 'image2D'
                        obj.TextureMap = in;
                     case {'uint8' 'double'}
                        obj.TextureMap = image2D('Image',in);
                     otherwise
                         error('Wrong input for Texture map');
                 end
                 %obj.TextureMap.Parent = obj;
                 checkVerticesColorUpdate(obj,'update texture');
        end
        function out = get.TextureMap(obj)
            out = obj.TextureMap;
            if~mySuperClass.isH(out), out = []; end
        end
        function obj = set.PoseLM(obj,poseLM)
                 if ~isempty(obj.PoseLM), delete(obj.PoseLM); end
                 if isempty(poseLM), obj.PoseLM = []; return; end
                 if ~strcmp(class(poseLM),'LMObj'), error('Input must be a LM object'); end
                 %if ~(size(poseLM.Location,2)==5), error('Exactly 5 pose landmarks are required'); end
                 poseLM.Parent = obj;
                 poseLM.Axes = obj.Axes;
                 obj.PoseLM = poseLM;
        end
        function out = get.PoseLM(obj)
            out = obj.PoseLM;
            if~mySuperClass.isH(out), out = []; end
        end
        function obj = set.CustomLM(obj,LM)
%                  if isempty(LM), obj.CustomLM = []; return; end
%                  if ~strcmp(class(LM),'LMObj'), error('Input must be a LM object'); end                 
%                  LM.Tag = 'Custom';
%                  LM.Parent = obj;
%                  obj.CustomLM = LM;
        end
        function obj = set.Area(obj,area)
                 % test whether Area class
%                  if ~strcmp(class(poseLM),'areaObj'), error('Input must be a area object'); end
%                  obj.Area = area;
%                  obj.Area.Parent = obj;
        end
        function A = get.Adjacency(obj)
            f = double(obj.Tri');
            A = sparse([f(:,1); f(:,1); f(:,2); f(:,2); f(:,3); f(:,3)], ...
                       [f(:,2); f(:,3); f(:,1); f(:,3); f(:,1); f(:,2)], ...
                                                                   1.0);
            % avoid double links
            A = double(A>0);
        end
        function nr = get.nrCLM(obj)
            nr = length(obj.CustomLM);
        end
        function nr = get.nrA(obj)
            nr = length(obj.Area);
        end
        function out = get.KDE(obj)
            out = obj.KDE;
            if ~mySuperClass.isH(out), out = []; return; end
        end
        function obj = set.KDE(obj,in)
                 if ~isempty(obj.KDE); delete(obj.KDE); end
                 obj.KDE = in; 
        end
        function obj = set.TextureColor(obj,Color)
                 if ~(size(Color,1)==3), Color = Color'; end
                 %if ~(size(Color,2)==size(obj.Vertices,2)), return; end
                 warning('off','All');
                 if max(Color(:))<=3
                    obj.TextureColor = uint8(Color*255);
                 else
                    obj.TextureColor = uint8(Color);
                 end
                 warning('on','All');
                 %if ~isempty(obj.TextureMap)&&~isempty(obj.UV), return; end
                 checkVerticesColorUpdate(obj,'update texture');
        end
        function obj = set.IndexedColor(obj,Color)
                 if ~(size(Color,1)==1), Color = Color'; end
                 obj.IndexedColor = Color;
                 checkVerticesColorUpdate(obj,'update index');
        end      
        function out = get.TextureColor(obj)
           if isempty(obj.TextureMap)||isempty(obj.UV), out = double(obj.TextureColor)/255; return; end
           out = getValues(obj.TextureMap,obj.UV);
           %out = obj.TextureColor;
        end
        function obj = set.UV(obj,in)
           obj.UV = in;
           checkVerticesColorUpdate(obj,'update texture');
        end
    end
    methods % InterFace functions
        function export(obj)
            [filename, pathname,filterindex] = uiputfile({ '*.mat','Mat Object';...
                                                           '*.mat','Mat Struct';...
                                                           '*.obj','Obj Wavefront';...
                                                           '*.wrl','Wrl';...
                                                           '*.ply', 'Ply';...
                                                           '*.smf', 'Smf';...
                                                           '*.off', 'Off';...
                                                            },'Export As',obj.Tag);
            if isequal([filename,pathname],[0,0]),return; end
            cd(pathname);
            switch filterindex
                case 1% Obj
                    exportObject(obj,filename(1:end-4));
                case 2% Obj
                    exportStruct(obj,filename(1:end-4));
                case 3% Obj
                    exportWavefront(obj,filename(1:end-4));
                      
                case 4% Wrl
                case 5% Ply
                case 6% Smf
                case 7% Off
                    exportOff(obj,filename(1:end-4));
            end            
        end
    end
    methods(Static = true)
        function obj = import(filename, varargin)
                 if nargin == 0
                    [nrfiles,filename,pathname,filterindex] = meshObj.importGetFile;
                    if isempty(filename), obj = []; return; end
                    obj = cell(1,nrfiles);
                    for i=1:1:nrfiles
                        obj{i} = meshObj;
                        obj{i}.Tag = filename{i}(1:end-4);
                        switch filterindex
                            case 'Mat Object'
                                importObject(obj{i},filename{i},'Path',pathname,'Type',filterindex);
                            case 'Mat Struct'
                                importStruct(obj{i},filename{i},'Path',pathname,'Type',filterindex);
                            case {'Wavefront' '3DMD 2pod' '3DMD 5pod' 'Eyetronics'}
                                importWavefront(obj{i},filename{i},'Path',pathname,'Type',filterindex);
                            case {'WavefrontMevisLab'}
                                importWavefrontMevisLab(obj{i},filename{i},'Path',pathname,'Type',filterindex);    
                            case {'FMM Original' 'FMM Mapped'}
                                importFMM(obj{i},filename{i},'Path',pathname,'Type',filterindex);
                            case 'Wrl'
                                importWrl(obj{i},filename{i},'Path',pathname,'Type',filterindex);
                            case 'Ply'
                                importPly(obj{i},filename{i},'Path',pathname,'Type',filterindex);
                            case 'Smf'
                                importSmf(obj{i},filename{i},'Path',pathname,'Type',filterindex);
                            case 'Off'
                                importOff(obj{i},filename{i},'Path',pathname,'Type',filterindex);
                            case 'Nrf'
                                importNrf(obj{i},filename{i},'Path',pathname,'Type',filterindex);
                            case 'Xls'
                                delete(obj{i});
                                obj{i}= LMObj;
                                obj{i}.Tag = filename{i}(1:end-4);
                                importXls(obj{i},filename{i},'Path',pathname,'Type',filterindex);
                                continue;
                            otherwise
                        end
                        border(obj{i});
                    end
                    if nrfiles == 1, obj = obj{1}; end% no need to return a cell array
                 else
                    switch filename(end-3:end)
                        case '.obj'
                            obj = meshObj;
                            obj.Tag = filename(1:end-4);
                            importWavefront(meshObj,filename,varargin{:});
                            border(obj);
                        case '.mat'
                        otherwise
                            error('Unknown file format');
                    end
                 end
        end
        function varargout = importGetFile
             [filename, pathname,filterindex] = uigetfile({'*.mat','Mat Object';...
                                                           '*.mat','Mat Struct';...
                                                           '*.obj','Obj Wavefront';...
                                                           '*.obj','Obj 3DMD(2Pod)';...
                                                           '*.obj','Obj 3DMD(5Pod)';...
                                                           '*.obj','Obj Eyetronics';...
                                                           '*.obj','Obj MevisLab';...
                                                           '*.mat', 'FMM original';...
                                                           '*.mat', 'FMM Mapped';...
                                                           '*.wrl','Wrl';...
                                                           '*.ply', 'Ply';...
                                                           '*.smf', 'Smf';...
                                                           '*.off', 'Off';...
                                                           '*.nrf', 'Nrf';...
                                                           '*.xls', 'Excell Landmarks';...
                                                            },'MultiSelect','on');
             if isequal([filename,pathname],[0,0]),varargout{1} = [];varargout{2}=[];varargout{3}=[];varargout{4} = [];return; end
             if ~iscell(filename)
                varargout{1} = 1;
                varargout{2} = {filename};
             else
                varargout{1} = size(filename,2);
                varargout{2} = filename;
             end
             varargout{3} = pathname;
             switch filterindex
                 case 1
                     varargout{4} = 'Mat Object';
                 case 2
                     varargout{4} = 'Mat Struct';
                 case 3
                     varargout{4} = 'Wavefront';
                 case 4
                     varargout{4} = '3DMD 2pod';
                 case 5
                     varargout{4} = '3DMD 5pod';
                 case 6
                     varargout{4} = 'Eyetronics';
                 case 7
                     varargout{4} = 'WavefrontMevisLab';    
                 case 8
                     varargout{4} = 'FMM Original';
                 case 9
                     varargout{4} = 'FMM Mapped';
                 case 10
                     varargout{4} = 'Wrl';
                 case 11
                     varargout{4} = 'Ply';
                 case 12
                     varargout{4} = 'Smf';
                 case 13
                     varargout{4} = 'Off';
                 case 14
                     varargout{4} = 'Nrf';
                 case 15
                     varargout{4} = 'Xls';
                 otherwise
             end
                %          cd(pathname);
        end
    end
end % classdef