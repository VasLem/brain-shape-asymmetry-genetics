function Mesh = fastrbf_normalsfrommesh( varargin )
% FASTRBF_NORMALSFROMMESH Create vertex normals from faces
%    M = FASTRBF_NORMALSFROMMESH(M) calculates vertex normals for
%    the mesh M based on the mesh faces.  The normals are placed in
%    the Gradient field of M.
%
%    M = FASTRBF_NORMALSFROMMESH(...,'attr','Name') returns the
%    output vertex normals in M.Name.Gradient instead of the default
%    M.Gradient.
%
%    See also: FASTRBF, FASTRBF_NORMALSFROMSIGNS, FASTRBF_NORMALSFROMSCAN,
%              FASTRBF_NORMALSFROMPOINTS

[alt,n] = FastRBF_MEX( 'NormalsFromMesh', varargin{:} );
Mesh = varargin{1};
if isempty(alt)
  Mesh.Gradient = n;
else
  Mesh = setfield(Mesh, alt, 'Gradient', n);
end
