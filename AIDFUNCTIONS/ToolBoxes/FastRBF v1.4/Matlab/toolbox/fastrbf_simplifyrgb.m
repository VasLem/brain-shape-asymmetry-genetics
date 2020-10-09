function M = fastrbf_simplifyrgb(varargin)
% FASTRBF_SIMPLIFYRGB Constrained simplification of an isosurface
%    M = FASTRBF_SIMPLIFYRGB(M,SOL,ERROR,R,G,B,RGBSTEP,ISOVALUE)
%    simplifies the isosurface mesh M so that it remains within 
%    ERROR distance of the ISOVALUE isosurface of the RBF solution SOL, 
%    and so that the RGB colour value does not change more than RGBSTEP 
%    along any edge.  RGB colour values are given in the three RBF 
%    solutions R, G and B.  The isosurface value of the surface is 
%    specified in ISOVALUE.  A new mesh is returned.
%    Specify -1 for ERROR to leave it unconstrained. One of 'faces' or
%    'points' must be specified in this case.
%
%    M = FASTRBF_SIMPLIFYRGB(M,SOL,ERROR,R,G,B,RGBSTEP) uses
%    ISOVALUE = 0.
%
%    FASTRBF_SIMPLIFYRGB(...,'faces',NF) will stop when the number
%    of faces drops to NF.  The default is 0, no limit.
%
%    FASTRBF_SIMPLIFYRGB(...,'points',NP) will stop when the number
%    of points drops to NP.  The default is 0, no limit.
%
%    FASTRBF_SIMPLIFYRGB(...,'angle',A) will set the minimum
%    enclosed angle to A.  The default is 10 degrees.
%
%    FASTRBF_SIMPLIFY(...,'inverted') specifies that the face orientation
%    has been inverted.
%
%    See also: FASTRBF, FASTRBF_SIMPLIFY

M = FastRBF_MEX('SimplifyRGB',varargin{:});
