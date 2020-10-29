function [X,Y,Z] = fastrbf_gridvec( Grid )
% FASTRBF_GRIDVEC Get coordinate vectors of grid values
%    [X,Y,Z] = FASTRBF_GRIDVEC(GRID) return vectors X, Y and Z 
%    such that the value at Grid[k,j,i] has the coordinates X[i], Y[j], Z[k].  
%
%    [X,Y] = FASTRBF_GRIDVEC(GRID) return vectors X and Y 
%    such that Grid[j,i] has the coordinates X[i], Y[j].
%
%    See also: FASTRBF, FastRBF_GridEval, FastRBF_GridCoords

siz = size(Grid.Value) - 1;
X = (0:siz(1)) * Grid.Spacing(1) + Grid.Min(1);
Y = (0:siz(2)) * Grid.Spacing(2) + Grid.Min(2);
if nargout == 3
  if length(siz) < 3
    Z = 0;
  else 
    Z = (0:siz(3)) * Grid.Spacing(3) + Grid.Min(3);
  end
end
