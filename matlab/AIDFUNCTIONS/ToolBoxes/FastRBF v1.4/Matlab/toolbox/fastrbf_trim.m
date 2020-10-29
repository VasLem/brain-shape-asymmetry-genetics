function p = fastrbf_trim( varargin )

% FASTRBF_Trim Trims scan, mesh or point list data to a set of points
%    Y = FASTRBF_TRIM(X, TRIMPTS, DIST) removes points from X that lie 
%    further than DIST from any point in TRIMPTS.  
%    Mesh faces that reference removed vertices are also removed.  
%    The Size field of scan objects are updated.
%
%    FASTRBF_TRIM(..., 'edges') trims only the edges of a mesh, 
%    thus leaving any holes filled in (unless they touch an edge).
%
%    See also: FASTRBF, FASTRBF_CROP

p = FastRBF_MEX( 'Trim', varargin{:} );
