function Points = fastrbf_normalsfrompoints( varargin )

% FASTRBF_NORMALSFROMPOINTS Create vertex normals from a point cloud
%    M = FASTRBF_NORMALSFROMPOINTS(P) estimates vertex normals for the
%    point list P. Note that in general, without additional information 
%    determining normals is ambiguous and using FASTRBF_NORMALSFROMSIGNS, 
%    FASTRBF_NORMALSFROMSCAN or FASTRBF_NORMALSFROMMESH is preferable.
%
%    M = FASTRBF_NORMALSFROMPOINTS(P, RADIUS)
%    RADIUS is the radius of the local patches for estimating
%    surface normals. It is a distance in the same units as P. A RADIUS of 0 
%    chooses the best plane from a series of patches where the size of each patch
%    is determined as the smallest radius which 
%    encloses 10, 30, 100 and 300 points. Default is 0.
%
%    M = FASTRBF_NORMALSFROMPOINTS(P, RADIUS, FACTOR) 
%    FACTOR is the minimum plane factor (or plane-y-ness) required for
%    the local patch for a normal to be accepted from it. A typical value
%    is 2, which means it's twice as much like a plane than a random 
%    cloud of points. FACTOR must be greater than 1. Default is 1.8.
%
%    M = FASTRBF_NORMALSFROMPOINTS(..., 'attr', 'Name') stores output
%    in M.Name.Gradient field instead of the default M.Gradient.
%    Input P may be a pointlist or a mesh.
%
%    See also: FASTRBF_NORMALSFROMSIGNS, FASTRBF_NORMALSFROMMESH,
%              FASTRBF_NORMALSFROMSCAN


[alt,n] = FastRBF_MEX( 'NormalsFromPoints', varargin{:} );
Points = varargin{1};
if isempty(alt)
  Points.Gradient = n;
else
  Points = setfield(Points, alt, 'Gradient', n);
end
