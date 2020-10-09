function [p,u] = fastrbf_unique( varargin )
% FASTRBF_UNIQUE Remove overlapping points
%    X = FASTRBF_UNIQUE(X,DIST,'absolute|relative') removes duplicate points from the
%    point list, mesh or scan X.  Points are considered identical
%    when the distance between them is less than DIST. Use 'absolute' 
%    or 'relative' to specify whether DIST is relative to the 
%    dimensions of the bounding box i.e. DIST*(longest side
%    of bounding box of points), or is an absolute spacing between points.
%
%    [X U] = FASTRBF_UNIQUE(X) returns a flag indicating if the
%     input points were unique.
%
%    X = FASTRBF_UNIQUE(X) uses DIST = 1e-10, as required by
%    FASTRBF_FIT.
%
%    See also: FASTRBF

if nargout == 2
  [p,u] = FastRBF_MEX( 'Unique', varargin{:} );
else
  p = FastRBF_MEX( 'Unique', varargin{:} );
end
