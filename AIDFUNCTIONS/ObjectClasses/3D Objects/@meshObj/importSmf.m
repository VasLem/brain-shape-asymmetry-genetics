function importSmf(obj,filename,varargin)
    % changing directory
         Input = find(strcmp(varargin,'Path'));
         if ~isempty(Input)
            cd(varargin{Input+1});
         end
     % Getting type
         Input = find(strcmp(varargin,'Type'));
         if isempty(Input)
            type = 'Off';
         else
            type = varargin{Input+1};
         end
     % quick inreading of vertices, normals and texturecoordinates;
         [vertex,face] = read_smf(filename);
      % Storing information in object
          obj.ColorMode = 'Single';
          obj.Vertices = vertex;
          obj.Faces = face;
      % storing kind of wavefront in Userdatas
          obj.UserData = lower(type);
end

function [vertex,face] = read_smf(filename)

% read_smf - read data from SMF file.
%
%   [vertex,face] = read_smf(filename);
%
%   'vertex' is a 'nb.vert x 3' array specifying the position of the vertices.
%   'face' is a 'nb.face x 3' array specifying the connectivity of the mesh.
%
%   Copyright (c) 2004 Gabriel Peyré

    fid = fopen(filename,'r');
    if( fid==-1 )
        error('Can''t open the file.');
        return;
    end

    vertex = [];
    face = [];

    str = 0;
    while ( str ~= -1)
        str = fgets(fid);   % -1 if eof
        if str(1)=='v'
            [a,str] = strtok(str);
            [a,str] = strtok(str); x = str2num(a);
            [a,str] = strtok(str); y = str2num(a);
            [a,str] = strtok(str); z = str2num(a);
            vertex = [vertex;[x y z]];
        elseif str(1)=='t' || str(1)=='f'
            [a,str] = strtok(str);
            [a,str] = strtok(str); x = str2num(a);
            [a,str] = strtok(str); y = str2num(a);
            [a,str] = strtok(str); z = str2num(a);
            face = [face;[x y z]];
        end
    end

    fclose(fid);
end