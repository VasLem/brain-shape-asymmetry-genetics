classdef BatchFile < superHandleClass
   properties
       FileName;
       MainFolder;
       SubFolder;
   end
   properties (Dependent = true)
       Extension;
       Name;
   end
   properties (Hidden = true, Dependent = true)% hidden properties but easy in the interaction
       FullFolder;
       FullName; 
   end
   methods %CONSTRUCTOR
       function obj = BatchFile(varargin)
          obj = obj@superHandleClass(varargin{:});
       end % myFile Constructor 
   end
   methods % GETTING
       function out = get.FullFolder(obj)
          out =  [obj.MainFolder obj.SubFolder];
       end
       function out = get.FullName(obj)
           out = [obj.FullFolder obj.FileName];
       end
       function out = get.Extension(obj)
          index = strfind(obj.FileName,'.');
          if isempty(index), out = []; return;end
          out = obj.FileName(index:end);
       end
       function out = get.Name(obj)
           index = strfind(obj.FileName,'.');
           if isempty(index), out = []; return;end
           out = obj.FileName(1:index-1);
       end
   end
   methods % SETTING
       function obj = set.FileName(obj,in)
           if contains(in,'.')
               obj.FileName = in;
           else
               obj.FileName = [in '.ukn'];% unknown format
           end
       end
       function obj = set.FullFolder(obj,in)
          in = BatchFile.addBackSlash(in);
          obj.MainFolder = in;
          obj.SubFolder = [];
       end
       function obj =set.MainFolder(obj,in)
           in = BatchFile.addBackSlash(in);
           obj.MainFolder = in;
       end
       function obj = set.SubFolder(obj,in)
           in = BatchFile.addBackSlash(in);
           obj.SubFolder = in;
       end
       function obj = set.Extension(obj,in)
           in = BatchFile.addFrontDot(in);
           obj.FileName = [obj.Name in];
       end
   end
   methods % INTERFACING
       function out = fileExist(obj)
           tmp = dir(obj.FullName);
           if isempty(tmp), out = false; else, out = true; end
       end
   end
   methods (Static = true)
       function out = addBackSlash(in)
                out = in;
                if isempty(out), return; end
                if ~strcmp(out(end),'/'), out = [out '/'];end
       end
       function out = addFrontDot(in)
                out = in;
                if isempty(out), return; end
                if~strcmp(out(1),'.'), out = ['.' out];end
       end
   end
end