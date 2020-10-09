function Sol = fastrbf_fit( varargin )
% FASTRBF_FIT Create an RBF solution that interpolates a density field
%    S = FASTRBF_FIT(P) returns an RBF solution that interpolates the 
%    P.Value values at the points in P.Location.  Per point fitting 
%    accuracy is specified with the 'Accuracy' field of P when performing a 'reduce'
%    fit. Per-point upper and lower error bars are specified with the 
%    'Lower' and 'Upper' fields of P when performing an errorbar fit.
%
%    Per-point accuracies, lower and upper errorbars are positive numbers representing
%    the absolute difference in the uncertainty of the value, e.g. for 9.8 +/- 0.1 the
%    per point accuracy is 0.1. For  9.8 +0.3, -0.2, the error bars are P.Lower=0.2 & P.Upper=0.3.  
%    When a per-point accuracy or lower or upper bound is present and a global accuracy is 
%    specified, the greater of the two at each point is taken.  
%    This allows the precision of a fit to be incrementally refined by controling the global 
%    accuracy. In the context of an errorbar fit the global accuracy is treated as a symmetric 
%    constraint on the value at each point. 
%   
%    FASTRBF_FIT(P,ACC) uses ACC as the default point accuracy.
%
%    FASTRBF_FIT(..., 'direct') uses the direct fitting algorithm.
%    All input points are used as centres. Best for data with little redundancy
%    relative to the fitting accuracy.
%    This is the default.
%
%    FASTRBF_FIT(..., 'reduce') uses the centre reduction fitting algorithm. 
%    Best for clean data with little or no discernable noise. 
%    This is the default.
%
%    FASTRBF_FIT(..., 'errorbar') Uses error bar fitting algorithm.  
%    If both the 'Lower' and the 'Upper' fields are present then the errorbar fitter will use 
%    asymetric bounds.  In all other cases the fit will use symmetric bounds.  If one of 
%    the 'Lower' or 'Upper' fields is present it will use it for the bounds.  If neither is 
%    present then the 'Accuracy' per-point field will be used if it is available.  
%    If a global accuracy (ACC) is present then it is used at every point as explained above.
%  
%    Best for noisy data or when the accuracy is quite loose.
%
%    FASTRBF_FIT(..., 'biharm') sets the basic function to the
%    bi-harmonic spline, which is the default.
%
%    FASTRBF_FIT(..., 'triharm') sets the basic function to
%    tri-harmonic.
%
%    FASTRBF_FIT(..., 'multi') sets the basic function to
%    generalized multiquadric. That is, phi(x) = (x^2+c^2)^{k/2}.
%
%    FASTRBF_FIT(..., 'c', CONST) sets the basic function constant to CONST.
%
%    FASTRBF_FIT(..., 'k', EXP) when used with 'multi', sets the 
%    multiquadric exponent. EXP must be odd and less than 6. Default is 1.
%
%    FASTRBF_FIT(..., 'degree', DEG) sets the RBF polynomial degree to DEG.
%    Use -1 for no polynomial. Minimum polynomial degree varies with basic
%    function. Default is 1, maximum is 2 (quadratic polynomial).
%
%    FASTRBF_FIT(..., 'rho', SM) sets the smoothing value to SM,
%    higher values will create a smoother fit.
%
%    FASTRBF_FIT(..., 'memlimit', m) Set a soft limit on the amount 
%    of RAM to use to m Mb.  May be exceeded if m is less than 
%    the core memory requirements.  More RAM means more speed.
%    The default is half the available physical RAM.
%
%    FASTRBF_FIT(..., 'confidence', C) sets the confidence in the
%    specified accuracy to C (0 < C <= 1).  That is, 100*C percent of
%    the points must have the specified accuracy.
%
%    FASTRBF_FIT(..., 'attr', 'Name') uses P.Name.Value instead of
%    P.Value for the point density values, and P.Name.Accuracy
%    instead of P.Accuracy for accuracy-per-point.
%
%    FASTRBF_FIT(..., 'initial', R) uses R as the initial guess
%    at the RBF solution.  
%
% The following options apply to the 'errorbar' fit type only. 
%
%    FASTRBF_FIT(..., 'noredundant') prevents the returned rbf from
%    containing redundant constraints.  The errorbar fitter solves
%    the convex quadratic programming problem: minimize some energy seminorm 
%    of interpolant, S, subject to f(i)-lower(i) <= S(i) <= f(i)+upper(i). 
%    The type of energy seminorm depends on the basic function.
%    "redundant constraint" means the lagrange multipliers for some 
%    constraints are negative. At the true solution of the convex QPP 
%    all the lagrange multipliers are positive.
%
%    FASTRBF_FIT(..., 'accfactor', ACCFACTOR) sets the 'exact' fitting 
%    accuracy as a fraction of the shortest error bar.  The solution 
%    to the above QPP is in fact an exact fit on a subset of the points
%    to the values f(i)-lower(i) or f(i)+upper(i).  Specifying a smaller
%    value will slow things down but return an rbf closer to the 
%    true minimum. (i.e. slightly less energy). The default is 0.1.
%    ACCFACTOR must be between 0.0 and 0.2.
%
%    FASTRBF_FIT(..., 'boundtypes', BT) specifies type of error bar
%    per point. 
%    If BT(k) = 1 the value at point k is constrained below.
%    If BT(k) = 2 the value at point k is constrained above.
%    If BT(k) = 3 the value at point k is constrained above and below.
%    Default is BT = [3,3,3,3,....]
%
%    See also: FASTRBF, FASTRBF_MAKERBF, FASTRBF_ENERGY

Sol = FastRBF_MEX( 'Fit', varargin{:} );
