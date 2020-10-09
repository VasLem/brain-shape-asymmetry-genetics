function importXls(obj,filename,varargin)
         path = readVarargin(varargin{:});
         cd(path);
         Vert = xlsread(filename);
         obj.Vertices = Vert';      
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