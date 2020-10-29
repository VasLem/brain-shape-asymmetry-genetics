function G = fastrbf_grideval( varargin )
% FASTRBF_GRIDEVAL Evaluate an RBF on a regular grid
%    G = FASTRBF_GRIDEVAL(S,...) evaluates the RBF solution S on the
%    grid specified in the remaining arguments.  The return value is a
%    FastRBF grid struct, with the RBF values in the 'Value' field.
%
%    G = FASTRBF_GRIDEVAL(S,ACC,...) sets the evaluation accuracy to
%    ACC.  The default accuracy is S.AchievedAcc/100. Must be present
%    if S.AchievedAcc is zero.
%
%    FASTRBF_GRIDEVAL(...,'min',MIN) sets the minimum corner of the
%    evaluation grid.
%
%    FASTRBF_GRIDEVAL(...,'max',MAX) sets the maximum corner of the
%    evaluation grid.
%
%    FASTRBF_GRIDEVAL(...,'spacing',SPACING) sets the spacing of the
%    grid.
%
%    FASTRBF_GRIDEVAL(...,'size',SIZE) sets the number of grid samples 
%    in each direction to SIZE.
%
%    Of the four options 'min', 'max', 'spacing' and 'size' three must
%    be present.  Specifying all four causes an error. 
%    If the grid is under specified then default values for
%    MIN and MAX are taken from S.DataMin and S.DataMax as required.
%
%    FASTRBF_GRIDEVAL(...,'gradient') also evaluates the RBF
%    gradient at each grid point.  The gradients are placed in the
%    Gradient field of the grid struct.
%
%    FASTRBF_GRIDEVAL(...,'smooth', WIDTH) filters the rbf values 
%    with a low pass filter of width WIDTH. 3D Only.
%
%    FASTRBF_GRIDEVAL(...,'pointlist') returns the RBF values in a
%    pointlist struct rather than a grid struct.
%
%    See also: FASTRBF, FASTRBF_GRIDCOORDS, FASTRBF_POINTEVAL

G = FastRBF_MEX( 'GridEval', varargin{:} );
