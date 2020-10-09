function X = fastrbf_load( varargin )
% FASTRBF_LOAD Load a FastRBF file
%    X = FASTRBF_LOAD(FILENAME) loads the file FILENAME and
%    returns the structure contained in that file.  The file may
%    contain a FastRBF point list, solution, grid, scan or mesh.
%    The file must be in ARANZ data file format.
%    
%    If a Location field is present the file is assumed to 
%    contain pointlist data and any fields ending in '_Value',
%    '_Accuracy', '_Lower', '_Upper' or '_Gradient' are created as 
%    substructs with the corresponding field name as the suffix.  
%    Essentially XXX_Value is replaced with XXX.Value.
%
%    See also: FASTRBF, FASTRBF_SAVE

X = FastRBF_MEX( 'Load', varargin{:} );
