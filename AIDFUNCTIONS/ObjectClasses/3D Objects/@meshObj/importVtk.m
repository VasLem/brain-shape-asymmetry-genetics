function importVtk(obj,filename,varargin)
    % changing directory
         Input = find(strcmp(varargin,'Path'));
         if ~isempty(Input)
            cd(varargin{Input+1});
         end
     % Getting type
         Input = find(strcmp(varargin,'Type'));
         if isempty(Input)
            type = 'Wrl';
         else
            type = varargin{Input+1};
         end
         [vertex,face] = read_vtk(filename,0);
      % Storing information in object
          %obj.ColorMode = 'Single';
          obj.Vertices = vertex;
          obj.Faces = face;
          %obj.ColorMode = 'Single';
      % storing kind of wavefront in Userdatas
          obj.UserData = lower(type);
end

function [vertex,face] = read_vtk(filename, verbose)

% read_vtk - read data from VTK file.
%
%   [vertex,face] = read_vtk(filename, verbose);
%
%   'vertex' is a 'nb.vert x 3' array specifying the position of the vertices.
%   'face' is a 'nb.face x 3' array specifying the connectivity of the mesh.
%
%   Copyright (c) Mario Richtsfeld

    if nargin<2
        verbose = 1;
    end

    fid = fopen(filename,'r');
    if( fid==-1 )
        error('Can''t open the file.');
        return;
    end

    str = fgets(fid);   % -1 if eof
    if ~strcmp(str(3:5), 'vtk')
        error('The file is not a valid VTK one.');    
    end

    %%% read header %%%
    str = fgets(fid);
    str = fgets(fid);
    str = fgets(fid);
    str = fgets(fid);
    nvert = sscanf(str,'%*s %d %*s', 1);

    % read vertices
    [A,cnt] = fscanf(fid,'%f %f %f', 3*nvert);
    if cnt~=3*nvert
        warning('Problem in reading vertices.');
    end
    A = reshape(A, 3, cnt/3);
    vertex = A;

    % read polygons
    str = fgets(fid);
    str = fgets(fid);

    info = sscanf(str,'%c %*s %*s', 1);

    if((info ~= 'P') && (info ~= 'V'))
        str = fgets(fid);    
        info = sscanf(str,'%c %*s %*s', 1);
    end

    if(info == 'P')

        nface = sscanf(str,'%*s %d %*s', 1);

        [A,cnt] = fscanf(fid,'%d %d %d %d\n', 4*nface);
        if cnt~=4*nface
            warning('Problem in reading faces.');
        end

        A = reshape(A, 4, cnt/4);
        face = A(2:4,:)+1;
    end

    if(info ~= 'P')
        face = 0;
    end

    % read vertex ids
    if(info == 'V')

        nv = sscanf(str,'%*s %d %*s', 1);

        [A,cnt] = fscanf(fid,'%d %d \n', 2*nv);
        if cnt~=2*nv
            warning('Problem in reading faces.');
        end

        A = reshape(A, 2, cnt/2);
        face = repmat(A(2,:)+1, 3, 1);
    end

    if((info ~= 'P') && (info ~= 'V'))
        face = 0;
    end

    fclose(fid);

    return
end