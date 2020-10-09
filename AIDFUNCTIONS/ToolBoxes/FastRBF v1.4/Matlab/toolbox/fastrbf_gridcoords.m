function [X,Y,Z] = fastrbf_gridcoords( Grid )
% FASTRBF_GRIDCOORDS Get coordinates of grid values.
%    [X,Y,Z] = FASTRBF_GRIDCOORDS(GRID) gets the x, y and z coordinates
%    of the elements in the grid structure such that GRID[i] is the value
%    corresponding to [X[i], Y[j], Z[k]].  If the grid is 2D then
%    the Z values will all be zero.
%
%    [X,Y] = FASTRBF_GRIDCOORDS(GRID) gets only the X and Y
%    coordinates of a 2D grid structure.
%
%    See also: FASTRBF

siz = size(Grid.Value);
if length(siz) == 2
  [X,Y] = ndgrid( range(Grid.Min(1), Grid.Spacing(1), siz(1)), ...
		  range(Grid.Min(2), Grid.Spacing(2), siz(2)) );
  if nargout == 3
    Z = zeros(size(X));
  end
elseif length(siz) == 3
  [X,Y,Z] = ndgrid( range(Grid.Min(1), Grid.Spacing(1), siz(1)), ...
		    range(Grid.Min(2), Grid.Spacing(2), siz(2)), ...
		    range(Grid.Min(3), Grid.Spacing(3), siz(3)) );
else
  error('GRID.Value should be 2D or 3D')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R = range(min,spacing,siz)

R = (0:(siz-1)) * spacing + min;
