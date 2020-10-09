classdef image2D < superClass
    properties
        Image = [];
        PixInterp = '*linear';
    end
    properties (Dependent = true)
        Dim;% Image dimension, third ; nr of Channels
        ImageSize;% X and Y size
        Res;
        Red; % Red Channel
        Green;% Green Channel
        Blue;% Blue Channel
        Gray;% Gray Channel
        uint8Image;
    end
    methods % Constructor
        function obj = image2D(varargin)
          obj = obj@superClass(varargin{:});
        end
    end
    methods % Special Setting & Getting
        function out = get.Dim(obj)
            out = size(obj.Image,3);
        end
        function out = get.Red(obj)
            if ~obj.Dim==3, out = []; return; end
            out = obj.Image(:,:,1);
        end
        function obj = set.Red(obj,in)
            obj.Image(:,:,1) = in;
        end
        function out = get.Green(obj)
            if ~obj.Dim==3, out = []; return; end
            out = obj.Image(:,:,2);
        end
        function obj = set.Green(obj,in)
            obj.Image(:,:,2) = in;
        end
        function out = get.Blue(obj)
            if ~obj.Dim==3, out = []; return; end
            out = obj.Image(:,:,3);
        end
        function obj = set.Blue(obj,in)
            obj.Image(:,:,3) = in;
        end
        function out = get.Gray(obj)
            if obj.Dim == 1, out = obj.Image; return; end
            out = rgb2gray(obj.Image); 
        end
        function obj = set.Image(obj,in)
            if strcmp(class(in),'uint8'), in = double(in)/255; end
            obj.Image = in;
            postImage(obj);
        end
        function out = get.ImageSize(obj)
            if isempty(obj.Image), out = []; return; end
            out = [size(obj.Image,2) size(obj.Image,1)];
        end
        function out = get.Res(obj)
            if isempty(obj.Image), out = []; return; end
            out = 1./(obj.ImageSize);
        end
        function out = get.uint8Image(obj)
            out = obj.Image;
            out = uint8(out*255);
        end
        function obj = set.uint8Image(obj,in)
            switch class(in)
                case 'uint8'
                    obj.Image = double(in)/255;
                otherwise
                    error('Input must be uint8');
            end
            
        end
        
    end
    methods % Interface functions
         function imshow(obj,Field)
             if nargin == 1
                 im = obj.Image;
             else
                 im = obj.(Field);
             end
             imshow(im);
         end
         function out = numel(obj)
             out = numel(obj.Image);
         end
         function postImage(obj) %#ok<INUSD>
             % dummy
         end
    end
    methods (Static = true)
        function obj = import(filename, varargin)
                 if nargin == 0
                    [nrfiles,filename,pathname] = image2D.importGetFile;
                    if isempty(filename), obj = []; return; end
                    obj = cell(1,nrfiles);
                    cd(pathname);
                    for i=1:1:nrfiles
                        obj{i} = image2D;
                        index = strfind(filename{i},'.');
                        obj{i}.Tag = filename{i}(1:index-1);
                        FMT = filename{i}(index+1:end);
                        obj{i}.Image = imread(filename{i},FMT);
                    end
                    if nrfiles == 1, obj = obj{1}; end% no need to return a cell array
                 else
                    obj = image2D;
                    index = strfind(filename,'.');
                    obj.Tag = filename(1:index-1);
                    FMT = filename(index+1:end);
                    obj.Image = imread(filename,FMT);    
                 end
        end
        function varargout = importGetFile
             [filename, pathname,filterindex] = uigetfile({'*.*','All Files';...
                                                           '*.bmp','BMP';...
                                                           '*.png','PNG';...
                                                           '*.jpg', 'Jpeg';...
                                                           '*.gif', 'gif';...
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
                %          cd(pathname);
        end
        
    end
end
% 
% classdef image2D < mySuperClass
%     properties
%         Image = [];
%         PixInterp = '*nearest';
%         %Tag = '';
% %         Gradient = [];
% %         Phase = [];
% %         Edge = [];
%     end
%     properties (Dependent = true)
%         Dim;% Image dimension, third ; nr of Channels
%         ImageSize;% X and Y size
%         Res;
%         Red; % Red Channel
%         Green;% Green Channel
%         Blue;% Blue Channel
%         Gray;% Gray Channel
%         uint8Image;
%     end
%     methods % Constructor
%         function obj = image2D(varargin)
%           if nargin > 0
%              Input = find(strcmp(varargin, 'Image'));if ~isempty(Input), obj.Image = varargin{Input+1}; end
% %              Input = find(strcmp(varargin, 'Gradient'));if ~isempty(Input), obj.Gradient = varargin{Input+1}; end
% %              Input = find(strcmp(varargin, 'Phase'));if ~isempty(Input), obj.Phase = varargin{Input+1}; end
% %              Input = find(strcmp(varargin, 'Edge'));if ~isempty(Input), obj.Edge = varargin{Input+1}; end
%           end
%         end
%     end
%     methods % Special Setting & Getting
%         function out = get.Dim(obj)
%             out = size(obj.Image,3);
%         end
%         function out = get.Red(obj)
%             if ~obj.Dim==3, out = []; return; end
%             out = obj.Image(:,:,1);
%         end
%         function obj = set.Red(obj,in)
%             obj.Image(:,:,1) = in;
%         end
%         function out = get.Green(obj)
%             if ~obj.Dim==3, out = []; return; end
%             out = obj.Image(:,:,2);
%         end
%         function obj = set.Green(obj,in)
%             obj.Image(:,:,2) = in;
%         end
%         function out = get.Blue(obj)
%             if ~obj.Dim==3, out = []; return; end
%             out = obj.Image(:,:,3);
%         end
%         function obj = set.Blue(obj,in)
%             obj.Image(:,:,3) = in;
%         end
%         function out = get.Gray(obj)
%             if obj.Dim == 1, out = obj.Image; return; end
%             out = rgb2gray(obj.Image); 
%         end
%         function obj = set.Image(obj,in)
%             if strcmp(class(in),'uint8'), in = double(in)/255; end
%             obj.Image = in;
%             postImage(obj);
%         end
%         function out = get.ImageSize(obj)
%             if isempty(obj.Image), out = []; return; end
%             out = [size(obj.Image,2) size(obj.Image,1)];
%         end
%         function out = get.Res(obj)
%             if isempty(obj.Image), out = []; return; end
%             out = 1./(obj.ImageSize);
%         end
%         function out = get.uint8Image(obj)
%             out = obj.Image;
%             out = uint8(out*255);
%         end
%         function obj = set.uint8Image(obj,in)
%             switch class(in)
%                 case 'uint8'
%                     obj.Image = double(in)/255;
%                 otherwise
%                     error('Input must be uint8');
%             end
%             
%         end
%         
%     end
%     methods % Interface functions
%          function copy(obj,cobj)
%             copy@mySuperClass(obj,cobj);
%             cobj.Image = obj.Image;
% %             cobj.Gradient = obj.Gradient;
% %             cobj.Phase = obj.Phase;
%          end
%          function struc = obj2struc(obj)
%              struc = obj2struc@mySuperClass(obj);
%              struc.Image = obj.Image;
% %              struc.Gradient = obj.Gradient;
% %              struc.Phase = obj.Phase;
%          end
%          function obj = struc2obj(obj,struc)
%              obj.Image = struc.Image;
% %              obj.Gradient = struc.Gradient;
% %              obj.Phase = struc.Phase;
%          end
%          function imshow(obj,Field)
%              if nargin == 1
%                  im = obj.Image;
%              else
%                  im = obj.(Field);
%              end
%              imshow(im);
%          end
%          function out = numel(obj)
%              out = numel(obj.Image);
%          end
%          function postImage(obj) %#ok<INUSD>
%              % dummy
%          end
%     end
%     methods (Static = true)
%         function obj = import(filename, varargin)
%                  if nargin == 0
%                     [nrfiles,filename,pathname] = image2D.importGetFile;
%                     if isempty(filename), obj = []; return; end
%                     obj = cell(1,nrfiles);
%                     cd(pathname);
%                     for i=1:1:nrfiles
%                         obj{i} = image2D;
%                         index = strfind(filename{i},'.');
%                         obj{i}.Tag = filename{i}(1:index-1);
%                         FMT = filename{i}(index+1:end);
%                         obj{i}.Image = imread(filename{i},FMT);
%                     end
%                     if nrfiles == 1, obj = obj{1}; end% no need to return a cell array
%                  else
%                     obj = image2D;
%                     index = strfind(filename,'.');
%                     obj.Tag = filename(1:index-1);
%                     FMT = filename(index+1:end);
%                     obj.Image = imread(filename,FMT);    
%                  end
%         end
%         function varargout = importGetFile
%              [filename, pathname,filterindex] = uigetfile({'*.*','All Files';...
%                                                            '*.bmp','BMP';...
%                                                            '*.png','PNG';...
%                                                            '*.jpg', 'Jpeg';...
%                                                            '*.gif', 'gif';...
%                                                             },'MultiSelect','on');
%              if isequal([filename,pathname],[0,0]),varargout{1} = [];varargout{2}=[];varargout{3}=[];varargout{4} = [];return; end
%              if ~iscell(filename)
%                 varargout{1} = 1;
%                 varargout{2} = {filename};
%              else
%                 varargout{1} = size(filename,2);
%                 varargout{2} = filename;
%              end
%              varargout{3} = pathname;
%                 %          cd(pathname);
%         end
%         
%     end
% end