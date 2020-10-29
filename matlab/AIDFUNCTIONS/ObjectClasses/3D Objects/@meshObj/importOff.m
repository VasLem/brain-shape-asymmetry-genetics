function importOff(obj,filename,varargin)
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
         [vertex,face] = read_off(filename);
      % Storing information in object
          obj.ColorMode = 'Single';
          obj.Vertices = vertex;
          obj.Faces = face;
      % storing kind of wavefront in Userdatas
          obj.UserData = lower(type);
end

function [vertex,face] = read_off(filename)

% read_off - read data from OFF file.
%
%   [vertex,face] = read_off(filename);
%
%   'vertex' is a 'nb.vert x 3' array specifying the position of the vertices.
%   'face' is a 'nb.face x 3' array specifying the connectivity of the mesh.
%
%   Copyright (c) 2003 Gabriel Peyré


    fid = fopen(filename,'r');
    if( fid==-1 )
        error('Can''t open the file.');
        return;
    end

    str = fgets(fid);   % -1 if eof
    if ~strcmp(str(1:3), 'OFF')
        error('The file is not a valid OFF one.');    
    end

    str = fgets(fid);
    [a,str] = strtok(str); nvert = str2num(a);
    [a,str] = strtok(str); nface = str2num(a);



    [A,cnt] = fscanf(fid,'%f %f %f', 3*nvert);
    if cnt~=3*nvert
        warning('Problem in reading vertices.');
    end
    A = reshape(A, 3, cnt/3);
    vertex = A;
    % read Face 1  1088 480 1022
    [A,cnt] = fscanf(fid,'%d %d %d %d\n', 4*nface);
    if cnt~=4*nface
        warning('Problem in reading faces.');
    end
    A = reshape(A, 4, cnt/4);
    face = A(2:4,:)+1;


    fclose(fid);

end