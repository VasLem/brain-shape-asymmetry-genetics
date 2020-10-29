function [X,RGB] = fastrbf_import( varargin )
% FASTRBF_IMPORT Import third-party files into FastRBF structures
%    X = FASTRBF_IMPORT(FILENAME) reads the file FILENAME,
%    determining the type of the file from the extension.  The return
%    value X could be a scan, mesh or 2D or 3D point list.
%
%    [X,RGB] = FASTRBF_IMPORT(...) returns rgb colour-per-vertex data
%    if it is present in the file.  Only OBJ format supports this.
%    If no colour data is present, RGB will be empty.
%
%    FASTRBF_IMPORT(...,'ptx') forces PTX format.
%    FASTRBF_IMPORT(...,'pto') forces PTO format.
%    FASTRBF_IMPORT(...,'obj') forces OBJ format.
%    FASTRBF_IMPORT(...,'stl') forces STL format.
%    FASTRBF_IMPORT(...,'txt') forces TXT format.
%
%    Formats:
%       PTX - Cyra scan format, returned as a scan.
%       PTO - Maptak/ARANZ scan format, returned as a scan.
%       OBJ - Alias|Wavefront format, returned as a 3D point list,
%             with Gradient values if 'vn' entries were found in the
%             file.  Colour-per-vertex data can be included via a
%             NONSTANDARD extension of writing red, green and blue
%             values (integer or float) after the point locations,
%             e.g.:
%                 v 5.768 2.586 1.002 0 128 255
%       STL - Stereolithography format (binary version only),
%             returned as a mesh structure.
%       TXT - Formatted text output, supports 2D or 3D point lists,
%             with all point list fields optional.  The 'format'
%             option controls how a line is read. The 'comment' option
%             tells fastRBF to skip lines which begin with the 
%             associated comment string   
%
%    See also: FASTRBF, FASTRBF_EXPORT

if nargout == 1
  X = FastRBF_MEX( 'Import', varargin{:} );
else
  [X, RGB] = FastRBF_MEX( 'Import', varargin{:} );
end
