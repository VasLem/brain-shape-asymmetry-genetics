function P = fastrbf_densityfromnormals( varargin )
% FASTRBF_DENSITYFROMNORMALS Create density data from a surface with normals
%    P = FASTRBF_DENSITYFROMNORMALS(Q,MIN,MAX) creates a density
%    field P from an input point list Q with normal vectors.  The
%    input points are all given density values of 0.  New points
%    are creating inside and outside the surface along the surface
%    normal vectors and have non-zero values.  The distance from the
%    surface of these new points is a least MIN and at most MAX.
%
%    FASTRBF_DENSITYFROMNORMALS(Q,MAX) uses MIN = 0.
%
%    FASTRBF_DENSITYFROMNORMALS(...,'onesided') will create new
%    points only on the outside of the surface.  The default is to
%    create points both inside and outside.
%
%    FASTRBF_DENSITYFROMNORMALS(...,'increase',INCR) will create at
%    most INCR*N new points, where N is the number of input points.
%    The default is 1, which doubles the number of points.  May be 0 if
%    'numgrid' is positive.
%
%    FASTRBF_DENSITYFROMNORMALS(..., 'fixedlen') uses a fixed projection
%    distance and changes the density value to be the distance to the
%    closest point, subject to certain conditions.  
%    MIN is ignored with this option.
%
%    FASTRBF_DENSITYFROMNORMALS(...,'numgrid', NUM) will add a grid
%    of approximately NUM off-surface points, subject to certain conditions.  
%    This helps prevent spurious surfaces occurring in the volume of interest.
%    The MIN and MAX arguments may be omitted if NUM is positive, implying
%    an 'increase' value of 0.  That is, only grid points are added.
%
%    FASTRBF_DENSITYFROMNORMALS(...,'attr', NAME) Get input from, and 
%    put output into the fields of attribute sub-struct called NAME 
%    rather than the point list.
%
%    See also: FASTRBF

P = FastRBF_MEX( 'DensityFromNormals', varargin{:} );
