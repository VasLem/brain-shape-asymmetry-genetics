function [OutMesh,H] = fastrbf_mesh2mesh(InMesh, Accuracy, Resolution, varargin)
% FASTRBF_MESH2MESH Use the FastRBF library for mesh repair
%    MESH = FASTRBF_MESH2MESH(MESH,ACC,RES) fits an RBF to an
%    existing (possibly partial) mesh, then isosurfaces that RBF to
%    produce a manifold and evenly sampled output mesh.
%    It gives warnings concerning errors in the input mesh but does not
%    correct these.
%
%    FASTRBF_MESH2MESH(...,'verbose') prints all messages.
%
%    FASTRBF_MESH2MESH(...,'quiet') prints only warning and error messages.
%
%    FASTRBF_MESH2MESH(...,'messages',N) sets the message level to N.
%
%    See also: FASTRBF
p = InMesh;

% Check mesh & print warning if an edge shares > 2 faces
fastrbf_checkmesh(p, varargin{:});

% find the minimum triangle edge length
p1 = p.Location(:,p.Tri(1,:));
p2 = p.Location(:,p.Tri(2,:));
p3 = p.Location(:,p.Tri(3,:));
distsq = [sum((p1-p2).^2, 1), sum((p2-p3).^2, 1), sum((p1-p3).^2, 1)];
meddist = sqrt(median(distsq));

% calculate normal lengths
minlength = Accuracy*2;
maxlength = max(meddist, Accuracy*10);

% generate normals from mesh faces
p = fastrbf_normalsFromMesh( p, varargin{:} );

% convert mesh with normals into points with density
p = fastrbf_densityFromNormals( p, minlength, maxlength, varargin{:} );

% fit
s = fastrbf_fit( p, Accuracy, varargin{:} );

% isosurf
o = fastrbf_isosurf( s, Resolution, varargin{:} );

% display outmesh
h = fastrbf_view( o, [] );

% output arguments
OutMesh = o;
if nargout == 2
  H = h;
end
