function M = fastrbf_simplify(varargin)
% FASTRBF_SIMPLIFY Constrained simplification of an isosurface mesh
%    MESH = FASTRBF_SIMPLIFY(M,SOL,ERROR,ISOVALUE) simplifies the
%    isosurface mesh M so that it remains within ERROR distance of the
%    ISOVALUE isosurface of the RBF solution SOL.  The isosurface value 
%    of the surface is specified in ISOVALUE.  A new mesh is returned.
%    Specify -1 for ERROR to leave it unconstrained. One of 'faces' or
%    'points' must be specified in this case.
%
%    FASTRBF_SIMPLIFY(M,SOL,ERROR) uses ISOVALUE = 0.
%
%    FASTRBF_SIMPLIFY(...,'faces',NF) will stop when the number of
%    faces drops to NF.  The default is 0, no limit.
%
%    FASTRBF_SIMPLIFY(...,'points',NP) will stop when the number of
%    points drops to NP.  The default is 0, no limit.
%
%    FASTRBF_SIMPLIFY(...,'angle',A) will set the minimum
%    enclosed angle in a simplified face to A.  The default is 10 degrees.
%
%    FASTRBF_SIMPLIFY(...,'inverted') specifies that the face orientation
%    has been inverted.
%
%    See also: FASTRBF, FASTRBF_SIMPLIFYRGB

M = FastRBF_MEX('Simplify',varargin{:});
