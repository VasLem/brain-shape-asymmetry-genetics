function v = fastrbf
% FASTRBF Fast 2D and 3D RBF toolbox.
%   
%   The FastRBF toolbox is a collection of RBF (radial basis
%   function) tools.  The full list of functions is given below.
%   Copyright (c) 1997-2001 FarField Technology Ltd
%
%
%Message Levels:
%   The FastRBF function produce controlled levels of messages.
%   The available message levels are:
%      1 - errors
%      2 - warnings
%      3 - success messages
%      4 - progress bar (default)
%      5 - informational messages
%
%   The mesage level can be set using the following options, which
%   are accepted by all FastRBF functions:
%      'quiet'          - set message level to 2
%      'verbose'        - set message level to 5
%      'messages', N    - set message level to N
%
%   For a description of the data structures used with this library
%   see FASTRBF_TYPES.
%
%See also:
%     Core:
%      · FASTRBF_FIT
%            Fit an rbf solution to a list of points and values.
%      · FASTRBF_POINTEVAL
%            Evaluate an rbf on a set of points.
%      · FASTRBF_GRIDEVAL
%            Evaluate an rbf on a regular grid.
%      · FASTRBF_ENERGY
%            Evaluate the energy in an RBF.
%
%     Data processing:
%      · FASTRBF_CROP
%            Remove data outside a given box.
%      · FASTRBF_GRIDCOORDS
%				 Returns X, Y and Z coordinates of a grid structure.
%      · FASTRBF_GRIDVEC
%				 Returns vectors X, Y and Z where XxYxZ are the grid points.
%      · FASTRBF_MAKERBF
%            Create a FastRBF solution struct from existing RBF.
%      · FASTRBF_MAKEPOINTLIST
%            Create a FastRBF point list struct from existing arrays.
%      · FASTRBF_TRANSFORMRBF
%				 Applys an affine transformation to an RBF (rigid body).
%      · FASTRBF_TRIM
%            Remove further than a given distance from a set of points.
%      · FASTRBF_UNIQUE
%            Remove overlapping points from a point list scan or mesh
%
%     Isosurface:
%      · FASTRBF_ISOSURF
%            Extract an iso-surface of an rbf solution
%      · FASTRBF_MCUBES
%            Extract an iso-surface from the values of a grid structure.
%      · FASTRBF_SIMPLIFY
%				 Simplify an iso-surface mesh.
%      · FASTRBF_SIMPLIFYRGB
%				 Simplify an iso-surface mesh preserving colour.
%     
%     Surface Fit:
%      · FASTRBF_CHECKMESH
%            Checks and fixes mesh face orientation consistency.
%      · FASTRBF_NORMALSFROMMESH
%            Calculate mesh vertex normals from face normals
%      · FASTRBF_NORMALSFROMSCAN
%            Estimate surface normals, starting from scan data
%      · FASTRBF_NORMALSFROMSIGNS
%            Estimate surface normals, starting from outward pointing vectors
%      · FASTRBF_NORMALSFROMPOINTS
%            Estimate normals from surface point data without any additional information
%      · FASTRBF_DENSITYFROMNORMALS
%            Convert points with normals into distance-to-surface
%            density data
%      · FASTRBF_MESH2MESH
%            A script to repair an existing mesh.
%
%     Visualisation:
%      · FASTRBF_VIEW
%            Display point list, mesh and scan objects
%
%     FastRBF Files:
%      · FASTRBF_SAVE
%      · FASTRBF_LOAD
%
%     Other File Formats:
%      · FASTRBF_EXPORT
%      · FASTRBF_IMPORT

% Copyright 2001 Applied Research Associates NZ Ltd
% Author: Tim Evans

if nargout == 0
	FastRBF_MEX('version');
else
	v = FastRBF_MEX('version');
end
