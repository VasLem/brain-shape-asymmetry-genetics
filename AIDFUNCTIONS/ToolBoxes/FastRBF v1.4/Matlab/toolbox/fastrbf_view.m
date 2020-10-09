function H = fastrbf_view( varargin )
% FASTRBF_VIEW Display FastRBF data structures
%    FASTRBF_VIEW(DATA,'type',COL) displays DATA using the display
%    type 'type' and colour COL.
%
%    DATA can be a FastRBF point list, mesh, scan, solution or grid.
%
%    'type' is a string that can contain the following characters:
%       'p' -> make point locations visible
%       'f' -> make faces visible, if faces are present in DATA
%       'e' -> make edge visible, if faces are present in DATA
%       'v' -> make the vectors DATA.Gradient visible as arrows
%       'n' -> uses DATA.Gradient as vertex normals, if 'f' is present
%
%    COL can be one of several different shapes (where Np is the
%    number of points):
%        1-by-Np -> indexed colour value per vertex
%        3-by-Np -> RGB colour per vertex
%        4-by-Np -> RGBA colour per vertex (Matlab 6 only)
%        3-by-1  -> RGB colour for patch
%        4-by-1  -> RGBA colour for patch (Matlab 6 only)
%        1-by-M  -> where M<N and M is whole specifies groups of
%                   points, each of which is given a separate
%                   indexed value
%
%    FASTRBF_VIEW(DATA,'type') uses a default value for COL.  This
%    default is:
%        if DATA.Size exists (i.e. DATA is a scan):
%            DATA.Size
%        if DATA.Value exists:
%            DATA.Value
%        else:
%            DATA.Location(end,:) (the height of the point)
%
%    FASTRBF_VIEW(DATA) uses a default value for 'type'.  If DATA is a
%    mesh or a 2D grid, the default is 'f', otherwise it is 'p'.
%
%    H = FASTRBF_VIEW(...) returns a vector of patch handles.  If the
%    type contains 'p', 'f', or 'e' then the first handle is the
%    point/face/edge patch.  If the type contains 'v' the remaining
%    handles will be the output of QUIVER or QUIVER3.
%
%    FASTRBF_VIEW(...,'--','param','value',...) specifies patch
%    param/value pairs for the main patch object.  The '--' argument
%    marks the start of these pairs.
%
%    See also: FASTRBF

% get arguments
[err,Data,Type,Colour,NDim,PArgs] = GetArgs(varargin);
error(err); % does nothing if err empty

% process arguments
if isempty(Type)
  Type = DefaultType(Data);
end

[err,Colour,alpha] = GetColour(Data,Colour);
error(err);

ax = newplot;

h = [];

if isempty(setdiff(Type,'n'))
  error('TYPE must contain at least on of ''pvfe''')
end

if ~isempty(intersect('pef',Type))
  params = GetParams(Data,Type,Colour,alpha);
  h(end+1) = patch( params, PArgs{:} );
end

if any('v' == Type)
  if ~isfield(Data,'Gradient')
    error('DATA.Gradient must exist for ''v'' view')
  end
  vh = VectorPatch(Data);
  h = [h;vh];
end

if nargout > 0
  H = h;
end

if ~ishold
  if NDim == 2 | NDim == 3
    view(NDim);
  end
  grid on;
  axis equal;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function H = VectorPatch(Data)

if size(Data.Location,1) == 2
  x = Data.Location(1,:);
  y = Data.Location(2,:);
  u = Data.Gradient(1,:);
  v = Data.Gradient(2,:);

  ih = ishold;
  if ~ih, hold on, end
  H = quiver(x,y,u,v);
  if ~ih, hold off, end
else
  x = Data.Location(1,:);
  y = Data.Location(2,:);
  z = Data.Location(3,:);
  u = Data.Gradient(1,:);
  v = Data.Gradient(2,:);
  w = Data.Gradient(3,:);

  ih = ishold;
  if ~ih, hold on, end
  H = quiver3(x,y,z,u,v,w);
  if ~ih, hold off, end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Params = GetParams(Data,Type,Colour,Alpha)

Np = size(Data.Location,2);
Nd = size(Data.Location,1);

have_alpha = ~isempty(Alpha);

% points
Params.Vertices = Data.Location';

% colour
Params.FaceVertexCData = Colour';
if have_alpha
  Params.FaceVertexAlphaData = Alpha';
end

% points
if any('p' == Type)
  Params.Marker = '.';
  Param.MarkerSize = 5;
end

% faces
if any('f' == Type) | any('e' == Type)
  % face-type patch
  if isfield(Data,'Tri')
    Params.Faces = Data.Tri;
  else
    Params.Faces = zeros(3,0);
  end
  if isfield(Data,'Quad')
     if ~isempty(findstr(version, 'R12.1')) & ~isempty(Params.Faces)
        % matlab 12.1 is broken wrt handling patches with different size faces
        % hence we triangulate using very simple algorithm stolen from reducepatch.m
        faces = Data.Quad';
        [frows fcols] = size(faces);
        if fcols > 3 
           newfaces = zeros(frows*(fcols-2),3);
           newfaces(:,1) = repmat(faces(:,1),fcols-2,1);
           for k = 2:fcols-1
              newfaces([1:frows]+frows*(k-2), 2:3) = [faces(:,k) faces(:,k+1)];
           end
           Params.Faces = [Params.Faces newfaces'];
        end
     elseif isempty(Params.Faces)
        Params.Faces = Data.Quad;
     else
        Params.Faces=[NaN+zeros(1,size(Params.Faces,2)); Params.Faces];
        nq = size(Data.Quad,2);
        Params.Faces(:,end+1:end+nq) = Data.Quad;
     end
  end
  Params.Faces = fliplr(Params.Faces');
else
  % point-type patch
  Params.Faces = [ 1:Np ; NaN+zeros(1,Np) ]';
end

% colours
if any('f' == Type)
  % show faces in colour-per-vertex
  Params.FaceColor = 'interp';
  if have_alpha
    Params.FaceAlpha = 'interp';
  end
  Params.EdgeColor = 'none';
  if any('e' == Type)
    % edges with faces: show in black
    Params.EdgeColor = 'black';
  end
  if any('p' == Type)
    % points with faces: show in black
    Params.MarkerEdgeColor = 'black';
    Params.MarkerFaceColor = 'black';
  end
else
  if any('p' == Type)
    % points without faces: make makers colour-per-vertex
    Params.MarkerEdgeColor = 'flat';
    Params.MarkerFaceColor = 'flat';
  end
  Params.EdgeColor = 'interp';
  if have_alpha
    Params.EdgeAlpha = 'interp';
  end
  if any('e' == Type)
    % edges without faces: show in interpolated colour-per-vertex
    Params.EdgeColor = 'interp';
    if have_alpha
      Params.EdgeAlpha = 'interp';
    end
    Params.FaceColor = 'white';
  else
    % no edges: hide edges
    Params.EdgeColor = 'none';
    Params.FaceColor = 'none';
  end
end

% vertex normals
if any('f' == Type) & any('n' == Type) & isfield(Data,'Gradient')
  Params.VertexNormals = Data.Gradient';
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Type = DefaultType(Data)

if isfield(Data,'Tri') | isfield(Data,'Quad')
  Type = 'f';
else
  Type = 'p';
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Error,Colour,Alpha] = GetColour(Data,Colour)

Alpha = [];
Error = '';

Np = size(Data.Location,2);

% default values if Colour is empty

if isempty(Colour)
%  if isfield(Data,'Size')
    % use scan-grouped colourx
%    Colour = Data.Size;
%  else
  if isfield(Data,'Value')
    % use Value field as colour
    Colour = Data.Value;
  else
    Colour = Data.Location(end,:);
  end
end

% change colour into valid FastVertexCData'

csiz = size(Colour);
if prod(csiz) == Np
  % valid index-per-vertex
  Colour = Colour(:)';
elseif csiz(1) == 4 & prod(csiz(2:end)) == Np
  % valid rgba-per-vertex
  Colour = reshape( Colour, 4, prod(csiz(2:end)) );
  Alpha = Colour(4,:);
  Colour = Colour(1:3,:);

  try
    get(0, 'DefaultSurfaceFaceAlpha');
    alpha_ok = 1;
  catch
    alpha_ok = 0;
  end
  if ~alpha_ok
    % can't do rgba
    Alpha = [];
    warning('rgba colour is not supported, ignoring')
  end
elseif csiz(1) == 3 & prod(csiz(2:end)) == Np
  % valid rgb-per-vertex
  Colour = reshape( Colour, 3, prod(csiz(2:end)) );
elseif csiz(1) == 3 & prod(csiz) == 3
  % valid rgb colour for whole patch
  Colour = reshape( Colour, 3, 1);
  Colour = repmat(Colour, 1, Np);
elseif csiz(1) == 4 & prod(csiz) == 4
  % valid rgba colour for whole patch
  Colour = reshape( Colour, 4, 1);
  Alpha = Colour(4);
  Colour = repmat(Colour(1:3), 1, Np);
  Alpha = repmat(Alpha, 1, Np);
  try
    get(0, 'DefaultSurfaceFaceAlpha');
    alpha_ok = 1;
  catch
    alpha_ok = 0;
  end
  if ~alpha_ok
    % can't do rgba
    Alpha = [];
    warning('rgba colour is not supported, ignoring')
  end
elseif length(csiz) == 2 & ...
       csiz(1) == 1 & csiz(2) < Np & ...
       (strncmp(class(Colour),'int',3) | all(round(Colour)==Colour))
  % grouped colour
  groups = Colour;
  n = prod(size(groups));
  Colour = n + 1 + zeros(1,Np);
  first = 1;
  for x=1:n
    last = first + groups(x) - 1;
    Colour(first:last) = x + zeros(1,groups(x));
    first = last + 1;
  end
else
  Error = 'invalid COLOUR argument';
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Error,Data,Type,Colour,NDim,PArgs] = GetArgs(Args)

% default values
Error = '';
Data = [];
Type = '';
Colour = [];
NDim = 0;

% split out patch arguments (all args following '--')
for cut=1:length(Args)
  if ischar(Args{cut}) & strcmp(Args{cut},'--')
    cut = cut - 1;
    break
  end
end
PArgs = Args(cut+2:end);
Args = Args(1:cut);
if ~isempty(PArgs) & bitand(length(PArgs), 1) ~= 0
  Error = 'extra arguments must be in ''param'',''value'' pairs';
  return
end

Error = nargchk(1,3,length(Args));
if ~isempty(Error)
  return
end

% get args
Data = Args{1};
if length(Args) == 2
  if isa(Args{2}, 'double')
    Colour = Args{2};
  else
    Type = Args{2};
  end
elseif length(Args) == 3
  Type = Args{2};
  Colour = Args{3};
end

% check args
if isstruct(Data)
  if isfield(Data,'Location')
    % valid point list
    NDim = size(Data.Location,1);
  elseif isfield(Data,'Centres') & isfield(Data,'Coeffs')
    % solution
    d.Location = Data.Centres;
    d.Value = Data.Coeffs(1:size(Data.Centres,2));
    Data = d;
    NDim = size(Data.Location,1);
  elseif isfield(Data,'Min') & isfield(Data,'Max') & isfield(Data,'Spacing') & isfield(Data,'Value')
    % grid result
    if ndims(Data.Value) == 2
      % 2D -> convert to points and faces
      Data = Grid2_Mesh(Data);
      NDim = 2;
    elseif ndims(Data.Value) == 3
      % 3D -> convert to points
      Data = Grid3_Points(Data);
      NDim = 3;
    end
  else
    Error = 'DATA should be a FastRBF point list, scan, mesh, solution or grid struct';
    return
  end
else
  [NDim,n] = size(Data);
  if ~any(NDim == [2 3])
    Error = 'DATA (as double array) should be 2-by-N or 3-by-N';
    return
  end
  Data = struct('Location', reshape(Data,NDim,n));
end

if ~isempty(Type) & ~ischar(Type)
  Error = 'TYPE must be a string';
  return
end

Type = unique(Type(:));
if ~isempty(setdiff(Type, 'efnpv'))
  Error = 'TYPE may contain only the characters ''efnpv''';
  return
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Mesh = Grid2_Mesh( Grid )

% a single face is [c, c+1, c+n, c+n+1]
% where:
%    c = a + n(b-1)
%    a in [1,n)
%    b in [1,m)

[x,y] = fastrbf_gridcoords(Grid);
Mesh.Location = [x(:) y(:) Grid.Value(:)]';
Mesh.Value = Grid.Value(:)';
[n,m] = size(Grid.Value);
[a,b] = ndgrid(1:(n-1), 1:(m-1));
a = a(:);
b = b(:);
c = a + n*(b-1);
Mesh.Quad = [ c, c+1, c+n+1, c+n ]';
if isfield(Grid,'Gradient')
  gsiz = size(Grid.Gradient);
  Mesh.Gradient = reshape( Grid.Gradient, gsiz(1), prod(gsiz(2:end)) );
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Points = Grid3_Points( Grid )

[x,y,z] = fastrbf_gridcoords(Grid);
Points.Location = [x(:) y(:) z(:)]';
Points.Value = Grid.Value(:)';
if isfield(Grid,'Gradient')
  gsiz = size(Grid.Gradient);
  Points.Gradient = reshape( Grid.Gradient, gsiz(1), prod(gsiz(2:end)) );
end

return
