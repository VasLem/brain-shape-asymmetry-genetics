function rbf = fastrbf_makerbf(centres, coeffs, basicfunc, c, k, polybase, degree)

% FASTRBF_MAKERBF Creates an RBF solution struct from existing data.
%    RBF = FASTRBF_MAKERBF(CENTRES, COEFFS) 
%    creates a FastRBF solution structure from the given data. 
%    CENTRES is DIM x N where DIM is 2 or 3.
%    COEFFS has N basic function coefficients followed by 
%    DIM+1 polynomial coefficients in the order 1, x, y, z 
%    for biharmonic spline basic function and (DIM+1)*(DIM+2)/2 
%    coefficients in the order 1, x, x^2, y, yx, y^2, z, zx, zy, z^2
%    for the tri-harmonic spline.
% 
%    FASTRBF_MAKERBF(CENTRES, COEFFS, 'biharm' | 'triharm' | 'multi') 
%    sets the basic function to the biharmonic, triharmonic and 
%    multiquadric basic functions respectively.
%
%    FASTRBF_MAKERBF(CENTRES, COEFFS, 'multi', c, k) sets the 
%    multiquadric parameters where phi(x) = (|x|^2 + c^2)^(k/2).
%    The default values are c = 0, k = 1.
%    
%    FASTRBF_MAKERBF(CENTRES, COEFFS, basicfunc, c, k, polybase, degree) 
%    sets the centre of the polynomial and it's degree. The default
%    polynomial centre is the origin. The degree is determined from the
%    number of polynomial coeffs in COEFFS.
%    
%    See also: FASTRBF, FASTRBF_FIT

if length(size(centres)) ~= 2
	error('centres must be a 2-dimensional array');
end

[dim n] = size(centres);
if dim ~= 2 & dim ~= 3
	if n == 2 | n == 3
		centres = centres';
	else
		error('centres must be a 2 x N or 3 x N array');
	end
end

% set defaults
if nargin < 3
	basicfunc = 0;
else
   if isequal(lower(basicfunc), 'biharm')
      basicfunc = 0;
   elseif isequal(lower(basicfunc), 'triharm')
	   basicfunc = 1;
   elseif isequal(lower(basicfunc), 'multi')
      basicfunc = 2;
   end
end

if nargin < 4
	c = 0;
end
if nargin < 5
   k = 1;
end
if nargin < 6
   polybase = zeros(1, dim);
end
if nargin < 7
   pdim = length(coeffs) - length(centres);
   if pdim == 0
      deg = -1;
   elseif pdim == 1
      deg = 0;
   elseif pdim == dim+1;
      deg = 1;
   elseif pdim == (dim+1)*(dim+2)/2
      deg = 2;
   else
      error('could not determine polynomial degree');
   end
end

rbf.AchievedAcc = 0;
rbf.DefaultEvalAcc = 1E-5;	% same as in c
rbf.Centres = centres;
rbf.Coeffs = coeffs(:)';
rbf.PolyBase = polybase;
rbf.DataMin = min(centres');
rbf.DataMax = max(centres');
rbf.BasicFunc = basicfunc;
rbf.BasicFuncParam = c;
rbf.BasicFuncParam2 = k;
rbf.PolyDegree = deg;
rbf.Rho = 0.0;
rbf.FitType = 0;
rbf.Version = 'FastRBF V1.4';
