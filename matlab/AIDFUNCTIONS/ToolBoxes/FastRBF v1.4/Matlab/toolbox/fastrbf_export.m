function fastrbf_export( varargin )
% FASTRBF_EXPORT Export FastRBF data to third-party file formats
%    FASTRBF_EXPORT(X, FILENAME) writes X to the file FILENAME
%    determining the type of the file from the extension.  The
%    input data X can be a mesh or a 2D or 3D point list, depending
%    on the file type.
%
%    FASTRBF_EXPORT(X,RGB,RGBMAX,FILENAME) specifies rgb values
%    for each vertex.  RGB is a 3-by-N array, and RGBMAX is a single
%    number specifying the maximum RGB value (usually either 1 or
%    255).  RGB colour is used, if present, for OBJ and VRML output.
%
%    FASTRBF_EXPORT(...,'dxf') forces DXF format.
%    FASTRBF_EXPORT(...,'obj') forces OBJ format.
%    FASTRBF_EXPORT(...,'stl') forces STL format.
%    FASTRBF_EXPORT(...,'txt') forces formatted text format.
%    FASTRBF_EXPORT(...,'vrml') forces VRML format.
%    FASTRBF_EXPORT(...,'iges') forces IGES format.
%
%    FASTRBF_EXPORT(...,'format',FMT) sets the formatted text
%    format sting to FMT.  The default is '%x %y %z\n'.  See the
%    full documentation for a description of this parameter.
%
%    FASTRBF_EXPORT(...,'precision',N) sets the number of decimal
%    places used in DXF output.  Default is full machine accuracy.
%
%    FASTRBF_EXPORT(...,'ascii') selects ASCII rather than binary
%    output for formats that support both (currently only STL).
%
%    FASTRBF_EXPORT(...,'nonurbs') selects type 106-12 for iges output.
%    Each face is written as a piecewise linear closed curve.
%
%    Formats:
%       DXF - AutoDesk DXF format, supports 3D point list or mesh.
%       OBJ - Alias|Wavefront format, supports meshes, with
%             optional surface normals vectors in the 'Gradient' field.
%             Colour-per-vertex data can be included via a NONSTANDARD
%             extension of writing red, green and blue values as
%             integers on the range 0..255 after the point locations,
%             e.g.:
%                v 5.768 2.586 1.002 0 128 255
%       STL - Stereolithography format, binary or ASCII, supports
%             meshes, but only saves triangular faces.
%       TXT - Formatted text output, supports 2D or 3D point lists,
%             with all point list fields optional.  The 'format'
%             option controls what is written.
%       VRML- VRML 2.0 scene file, supports meshes.  The mesh is
%             written as a single object, with automatically
%             determined camera position.  RGB colour is supported.
%       IGES- IGES 5.1 (v9) graphic exchange format, supports points (type 116)
%             and meshes (types 128 and 106-12). By default each face
%             is written as a (flat) nurb patch (type 128).
%
%    See also: FASTRBF, FASTRBF_IMPORT

FastRBF_MEX( 'Export', varargin{:} );
