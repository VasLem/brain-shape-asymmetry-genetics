function fastrbf_save( varargin )
% FASTRBF_SAVE Save a FastRBF structure to a FastRBF file
%    FASTRBF_SAVE(X,FILENAME) saves the structure X in the file
%    FILENAME.  X can be a point list, solution, grid, mesh or scan.
%    The file is saved using the ARANZ data file format.
%
%    FASTRBF_SAVE(...,'ascii') saves in ASCII format rather than
%    the default binary.
%
%    See also: FASTRBF, FASTRBF_LOAD

FastRBF_MEX( 'Save', varargin{:} );
