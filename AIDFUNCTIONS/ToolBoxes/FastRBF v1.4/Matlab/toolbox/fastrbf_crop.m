function p = fastrbf_crop( varargin )

% FASTRBF_CROP Crops scan, mesh or point list data
%    Y = FASTRBF_CROP(X, MIN, MAX) removes points from X that lie 
%    outside the region specified by MIN and MAX.  Mesh faces that 
%    reference removed vertices are also removed.  
%    The Size field of scan objects are updated.
%
%    Y = FASTRBF_CROP(SCAN, MIN, MAX, 'range') removes points from
%    SCAN based on distance to the appropriate scan origin.  
%    Points closer than MIN and further than MAX are removed.
%
%    See also: FASTRBF, FASTRBF_TRIM

p = FastRBF_MEX( 'Crop', varargin{:} );
