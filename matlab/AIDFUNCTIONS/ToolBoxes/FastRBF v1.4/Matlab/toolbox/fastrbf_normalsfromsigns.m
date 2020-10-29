function Points = fastrbf_normalsfromsigns( varargin )

% FASTRBF_NORMALSFROMSIGNS Create vertex normals from scan
%    M = FASTRBF_NORMALSFROMSIGNS(P) estimates vertex
%    normals for the point list P, using outward pointing (but not
%    normal) vectors for initial data.  The outward pointing vectors
%    should be in the Gradient field of P and they will be replaced by
%    the actual normal vectors, or zero vectors where a normal could
%    not be reliably determined.
%
%    M = FASTRBF_NORMALSFROMSIGNS(P, RADIUS)
%    RADIUS is the radius of the local patches for estimating
%    surface normals. May be 0 to use the best plane from patches
%    with 10, 30, 100 and 300 points. Default is 0.
%
%    M = FASTRBF_NORMALSFROMSIGNS(P, RADIUS, FACTOR) 
%    FACTOR is the minimum plane factor (or plane-y-ness) required for
%    the local patch for a normal to be accepted from it. A typical value
%    is 2, which means it's twice as much like a plane than a random 
%    cloud of points. FACTOR must be greater than 1. Default is 1.8.
%
%    M = FASTRBF_NORMALSFROMSIGNS(..., 'attr', 'Name') uses the
%    vertex signs in M.Name.Gradient field instead of the default
%    M.Gradient.
%
%    See also: FASTRBF_NORMALSFROMMESH, FASTRBF_NORMALSFROMPOINTS,
%              FASTRBF_NORMALSFROMSCAN


[alt,n] = FastRBF_MEX( 'NormalsFromSigns', varargin{:} );
Points = varargin{1};
if isempty(alt)
  Points.Gradient = n;
else
  Points = setfield(Points, alt, 'Gradient', n);
end
