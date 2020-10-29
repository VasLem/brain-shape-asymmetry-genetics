function MeshInfo = fastrbf_checkmesh( varargin )
% FASTRBF_CHECKMESH Check mesh has consistenly ordered faces
%    INFO = FASTRBF_CHECKMESH(M) checks the mesh has consistently
%    ordered faces (based on edge-vertex connectivity) and is manifold*.
% 
%    INFO = FASTRBF_CHECKMESH(M, PERPART) calculates Consistent, 
%    Manifold and Closed properties for each part if PERPART is non-zero.
%
%    INFO is a struct with fields NumParts, Consistent, Manifold, Closed
%    and PartSizes. The last 3 values are boolean. PartSizes is a 
%    4-by-NumParts array, PartSizes(1,:) contains the assigned part numbers,
%    PartSizes(2,:) contains the number of faces per part.
%    PartSizes(3,:) contains the number of triangles per part.
%    PartSizes(4,:) contains the number of quadrilaterals per part.
%    If a mesh consists of multiple parts then all parts must be closed 
%    for the Closed field to be true.
%
%    See also: FASTRBF

MeshInfo = FastRBF_MEX( 'CheckMesh', varargin{:} );

% replace bit-coded properties with fields.
pp = MeshInfo.Properties;
MeshInfo.Consistent = (bitand(pp, 1) > 0);
MeshInfo.Manifold   = (bitand(pp, 2) > 0);
MeshInfo.Closed     = (bitand(pp, 4) > 0);
MeshInfo = rmfield(MeshInfo, 'Properties');

% replace per-part bit-coded properties with fields.
if isfield(MeshInfo, 'PartProperties')
   pp = MeshInfo.PartProperties;
	PartProp.Consistent = (bitand(pp, 1) > 0);
	PartProp.Manifold   = (bitand(pp, 2) > 0);
   PartProp.Closed     = (bitand(pp, 4) > 0);
   MeshInfo.PartProperties = PartProp;
end

