function Scan = fastrbf_normalsfromscan( varargin )

% FASTRBF_NORMALSFROMSCAN Create vertex normals from scan
%    M = FASTRBF_NORMALSFROMSCAN(S) estimates vertex
%    normals for the scan S, using the scanner location data to
%    provide initial information.  The normals are placed in the
%    Gradient field of M.
%
%    M = FASTRBF_NORMALSFROMSCAN(S, RADIUS)
%    RADIUS is the radius of the local patches used for estimating
%    surface normals.  May be 0 to use the best plane from 
%    10, 30, 100 or 300 points per patch. Default is 0.
%
%    M = FASTRBF_NORMALSFROMSCAN(S, RADIUS, FACTOR)
%    FACTOR is the minimum plane factor (or plane-y-ness) required for
%    the local patch for a normal to be accepted from it. A typical value
%    is 2, which means it's twice as much like a plane than a random 
%    cloud of points. FACTOR must be greater than 1. Default is 1.8.
%
%    M = FASTRBF_NORMALSFROMSCAN(..., 'attr', 'Name') returns the
%    output vertex normals in M.Name.Gradient instead of the default
%    M.Gradient.
%
%    See also: FASTRBF_NORMALSFROMSIGNS, FASTRBF_NORMALSFROMPOINTS,
%              FASTRBF_NORMALSFROMMESH


[alt, n] = FastRBF_MEX( 'NormalsFromScan', varargin{:} );
Scan = varargin{1};
if isempty(alt)
  Scan.Gradient = n;
else
  Scan = setfield(Scan, alt, 'Gradient', n);
end
