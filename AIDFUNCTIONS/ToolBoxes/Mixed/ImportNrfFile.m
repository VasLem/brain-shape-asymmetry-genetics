function [mesh] = ImportNrfFile(filename)


[x, y, z, RGB, v, f, c] = nrfShow2(filename);

mesh.Location(1,:) = x;
mesh.Location(2,:) = y;
mesh.Location(3,:) = z;

if size(f,2) == 4
    mesh.Quad = f';
else
    mesh.Tri = f';
end

mesh.Color.Value = c';
end


function [x, y, z, RGB, v, f, c] = nrfShow2(fName)
% Opens an NEC .nrf (Fiore 3D digitiser) image file and parses it.
% First reads headerSize bytes (we know that a .nrf file starts with
% 'NECRF NULL' + 16 bytes which seem to be four 32-bit integers: 1 22 640 640
% The rest of the file consists of x,y,z triples of IEEE 32 bit
% floats, each followed by B,G,R,0 (zero) in four bytes.
% nrfShow Returns x, y and z arrays of data, an associated colour image and the
% vertex, faces and colour arrays
% displays the file as a Matlab figure

% if exist(fName)
    fileID = fopen(fName,'r');
% else
    %errordlg(['File ' fName ' does not exist']);
%     return;
% end
headerSize = 6; 
headerBytes = fread(fileID, headerSize, 'uchar');% 'NECRF...'
sixteenbytes = fread(fileID, 4, 'int32')'; % read the next 16 bytes
% seems always to be 1    22   640   640 if read as 32-bit integers

xyz = zeros(640,640,3);
texture = zeros(640,640,3);
count = 1;
for row = 1:640
    for column = 1:640 
        xyz(row, column, :) = fread(fileID, 3, 'float32')'; % x,y,z (note the transpose on the output of the fread
        texture(row, column, :) = fread(fileID, 3, 'uchar')'; % R,B,G (note the transpose...)
        spare = fread(fileID, 1, 'uchar'); % extra colour byte (alpha channel?, probably not used)
    end
end

% process RGB to suit Matlab plotting functions
RGB(:,:,1) = texture(:,:,3);  % interchange colours to suit MatLab R,G,B order
RGB(:,:,2) = texture(:,:,2);
RGB(:,:,3) = texture(:,:,1);
RGB = RGB ./ 256; % scale for Matlab imshow
c = reshape(RGB,(640*640),3); % process RGB to suit Matlab patch function

% process xyz to suit Matlab patch function
xyz(xyz==-100000) = NaN; % NEC use x,y,z = -100000 for data points off the surface of the face
x=xyz(:,:,1);
y=xyz(:,:,2);
z=xyz(:,:,3);
x=reshape(x,1,(640*640));
y=reshape(y,1,(640*640));
z=reshape(z,1,(640*640));

% process coordinates to suit "patch" v is the array of vertices
v(:,1)=x' .* -1; % flip left for right so rt side of face is -ve x
v(:,2)=y' .* -1; % flip up for down so that the top of the head is +ve y
v(:,3)=z';

f = makeFaces(640); % make the array of face connectivity info for "patch"

% display the 3D data as a patch
%figure;
%fvc.vertices = v;
%fvc.faces = f;
%fvc.facevertexcdata = c;
%fvcr = reducepatch(fvc, 0.01);
%patchHandle = patch(fvc);
%figure;
%patch(fvcr);
% patchHandle = patch('Vertices',v,'Faces',f,'FaceVertexCData',c);
% axis equal;
% axis tight;
% 
% axesHandle = gca;
% set(axesHandle, 'CameraPosition', [2553 475 -2127],...
%     'CameraTarget', [0 0 0],...
%     'CameraUpVector', [0 1 0]); % set +ve y as up direction
% 
% set(patchHandle, 'BackFaceLighting', 'unlit',...
%     'EdgeColor', 'interp',...
%     'EdgeLighting', 'none',...
%     'FaceColor', 'interp',...
%     'FaceLighting', 'flat');
% 
fclose(fileID);
end

function faceData = makeFaces(n)
% generates a matrix of vertex connectivity info for an n x n matrix of
% vertex coordinates - faceData. faceData is an (n-1)^2 x 4 list of vertex
% numbers. Each group of 4 representing one facet.

nMinus1 = n-1;
faceData = zeros(nMinus1*nMinus1, 4);
for rowNo = 1:n-1
    for colNo = 1:n-1
        i = colNo + rowNo*n - n;
        faceData((nMinus1 * (rowNo - 1) + colNo),:) = [i (i+1) (i+1+n) (i+n)];
    end
end
end