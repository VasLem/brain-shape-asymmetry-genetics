function Mesh = fastrbf_mcubes( varargin )
% FASTRBF_MCUBES Extract an isosurface from a 3D grid
%    M = FASTRBF_MCUBES(GRID,ISOVALUE) extracts an isosurface
%    with isosurface value ISOVALUE from the values in GRID
%    The isosurface is returned in the FastRBF mesh object M.
%
%    FASTRBF_MCUBES(GRID) uses ISOVALUE = 0.
%
%    FASTRBF_MCUBES(...,'invert') inverts face orientation.
%
%    FASTRBF_MCUBES(...,'open'|'closeminus'|'closeplus') specifies
%    boundary capping. 'open', the defualt, does not cap the 
%    iso-surface at the boundaries.  'closeminus' caps at boundary
%    and assigns points outside the box a value less than the 
%    iso-surface value. 'closeplus' caps at the boundary
%    and assigns points outside the box a value above than the 
%    iso-surface value.
%
%    FASTRBF_MCUBES(...,'tri'|'triquad') sets the type of
%    faces generated: triangles, or both triangles and quadrilaterals.
%    The default is 'tri'.
%
%    See also: FASTRBF, FASTRBF_VIEW


Mesh = FastRBF_MEX( 'MCubes', varargin{:} );
