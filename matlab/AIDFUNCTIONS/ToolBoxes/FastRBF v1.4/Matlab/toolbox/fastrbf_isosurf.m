function Mesh = fastrbf_isosurf( varargin )
% FASTRBF_ISOSURF Extract an isosurface from an RBF solution
%    M = FASTRBF_ISOSURF(S,RES,ISOVALUE,ACC) extracts an isosurface
%    with isosurface value ISOVALUE and resolution RES from the RBF
%    solution S.  The RBF evaluation accuracy is ACC.  The isosurface
%    is returned in the FastRBF mesh object M.  The default bounding
%    box of the isosurface is the bounding box of the input RBF S.
%
%    M = FASTRBF_ISOSURF(S,RES,ISOVALUE) uses ACC = S.DefaultEvalAcc
%
%    M = FASTRBF_ISOSURF(S,RES) uses ISOVALUE = 0
%
%    FASTRBF_ISOSURF(...,'min',MIN) sets the minimum corner of the
%    evaluation box.
%
%    FASTRBF_ISOSURF(...,'max',MAX) sets the maximum corner of the
%    evaluation box.
%
%    FASTRBF_ISOSURF(...,'fit',FIT) tells the isosurfacer to 
%    adjust the sampling grid to fit exactly to the bounding box.
%    The characters 'x', 'y' and 'z' in FIT specify
%    which axes to adjust along.  The default is FIT = '', which
%    uses a regular tetrahedral grid.
%
%    FASTRBF_ISOSURF(...,'normals') enables generation of
%    approximate vertex normals.
%
%    FASTRBF_ISOSURF(...,'invert') inverts face orientation.
%
%    FASTRBF_ISOSURF(...,'open'|'closeminus'|'closeplus') specifies
%    boundary capping. 'open', the defualt, does not cap the 
%    iso-surface at the boundaries.  'closeminus' caps at boundary
%    and assigns points outside the box a value less than the 
%    iso-surface value. 'closeplus' caps at the boundary
%    and assigns points outside the box a value above than the 
%    iso-surface value.
%
%    FASTRBF_ISOSURF(...,'optimal'|'plane'|'grid') sets the face
%    layout.  'optimal' uses full optimisation, 'plane' will not move
%    face edges off x,y-planes, and 'grid' will not move any edges,
%    preserving the tetrahedral grid.  The default is 'optimal'.
%
%    FASTRBF_ISOSURF(...,'orient', T) transforms the evaluation
%    lattice by the 4x4 matrix T.  This can be used with 'plane' 
%	  to set the orientation of the parallel planes. 
%
%    FASTRBF_ISOSURF(...,'up', DIR) simplified version of 'orient'.
%    DIR is the normal vector of the plane.  DIR does not need to 
%    have unit length.
%
%    FASTRBF_ISOSURF(...,'bboxfixed') prevents the isosurfacer
%    from transforming the bounding box when 'orient' or 'up' is used.
%
%    FASTRBF_ISOSURF(...,'tri'|'quad'|'triquad') sets the type of
%    faces generated: triangles, quadrilaterals or both.  The
%    default is 'tri'.
%
%    FASTRBF_ISOSURF(...,'seeds',P) sets the surface following
%    seeds to the point list locations in P.  The default is to use
%    the RBF centres.
%
%    FASTRBF_ISOSURF(...,'seedlimit',N) sets the number of seeds to
%    use to N.  The default is 1.  Higher numbers are slower, but are
%    required in order to find other, separate parts of the surface.
%    A negative or zero value will use all available seeds.
%
%    FASTRBF_ISOSURF(...,'smooth', WIDTH) filters the rbf values 
%    with a low pass filter of width WIDTH.
%
%    See also: FASTRBF, FASTRBF_SIMPLIFY, FASTRBF_VIEW


Mesh = FastRBF_MEX( 'Isosurf', varargin{:} );
