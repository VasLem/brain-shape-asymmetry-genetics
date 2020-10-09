function [a,b] = fbracket(funfcn,x,options,varargin)
%FBRACKET Scalar nonlinear function minimum bracketing
%   [X1,X2] = FBRACKET(FUN,x) returns values X1 and X2 that bound the local minimum near X
%   of the function that is described in FUN (usually an M-file). The function FUN should
%   return a scalar function value F when called with feval: F=feval(FUN,X).  See below
%   for more options for FUN. The values X1 and X2 can then be used as input to FMINBND.
%
%   [X1,X2] = FBRACKET(FUN,x,OPTIONS) minimizes with the default optimization
%   parameters replaced by values in the structure OPTIONS, created with
%   the OPTIMSET function.  See OPTIMSET for details.  FBRACKET uses these options: 
%   BracketStep, MaxBracketStep, MultipleBracketStep, MultipleBracketStepRange
%   Do 'type fbracket' for help...
%
%   [X1,X2] = FBRACKET(FUN,x,OPTIONS,P1,P2,...) provides for additional
%   arguments, which are passed to the objective function, FUN(X,P1,P2,...).
%   (Use OPTIONS = [] as a place holder if no options are set.)
%
%   See also OPTIMSET, FMINBND. 

%   Reference: "Numerical Recipes in C, 2nd edition"
%   Press et al, 1992

%   Original coding by Frederik Maes, Katholieke Universiteit Leuven, November 2000


% PARAMS
% Default values
step = 1;  
maxstep = 10;
multistep = 0;
multisteprange = [-10 -5 -1 1 5 10];
% These can be changed by specifying the proper options
if nargin>2 && ~isempty(options)
  if isfield(options,'BracketStep')
    step = options.BracketStep;    
  end
  if isfield(options,'MaxBracketStep')
    maxstep = options.MaxBracketStep;
  end
  if isfield(options,'MultipleBrackStep')
    multistep = options.MultipleBracketStep;
  end
  if isfield(options,'MultipleBracketStepRange')
    multisteprange = options.MultipleBracketStepRange;
  end
end
% CONSTANTS
MI_GOLD = 1.618033989;
MI_TINY = eps;

% Convert to inline function as needed
funfcn = fcnchk(funfcn,length(varargin));

% Initialize bracket (a,b,c) such that f(a)>f(b)
% First point is current point
ax = 0;
fa = feval(funfcn,ax,varargin{:});
% Second point
if ~multistep
  % Single point
  bx = ax + step;
  fb = feval(funfcn,bx,varargin{:});
elseif multistep
  % Test several second points
  bx = ax + step.*multisteprange;
  fb = [];
  for k=1:length(bx)
    fb(k) = feval(funfcn,bx(k),varargin{:});
  end
  % Select 'best' one
  [fb,k] = min(fb);
  bx = bx(k);
end
% Order 
if fb>fa
  dum = ax; ax = bx; bx = dum;
  dum = fa; fa = fb; fb = dum;
end
% Third point
cx = bx + MI_GOLD * (bx-ax);
fc = feval(funfcn,cx,varargin{:});


% Update bracket (a,b,c) such that f(a)>f(b) and f(c)>f(b)
while (fb > fc)
  r = (bx-ax) * (fb-fc);
  q = (bx-cx) * (fb-fa);
  t = sign(max(abs(q-r),MI_TINY)*(q-r));
  u = bx - ((bx-cx)*q-(bx-ax)*r) / (2.0*t*max(abs(q-r),MI_TINY));
  ulim = bx + maxstep * (cx-bx);
  fu = 0;
  if ((bx-u)*(u-cx) > 0.0)  
    fu = feval(funfcn,u,varargin{:});
    if (fu < fc)
      ax = bx; bx = u;
      fa = fb; fb = fu;
      break
    elseif (fu > fb)
      cx = u; 
      fc = fu;
      break
    end      
    u = cx + MI_GOLD * (cx-bx);
    fu = feval(funfcn,u,varargin{:});
  elseif ((cx-u) * (u-ulim) > 0.0)
    fu = feval(funfcn,u,varargin{:});
    if (fu < fc)
      bx = cx; cx = u;
      fb = fc; fc = fu;
      u = cx + MI_GOLD * (cx-bx);
      fu = feval(funfcn,u,varargin{:});
    end
  elseif ((u-ulim)*(ulim-cx) >= 0.0)
    u = ulim;
    fu = feval(funfcn,u,varargin{:});
  else
    u = cx + MI_GOLD * (cx-bx);
    fu = feval(funfcn,u,varargin{:});
  end
  
  ax = bx; bx = cx; cx = u;
  fa = fb; fb = fc; fc = fu;
end

% Reorder 
a = min(ax,cx);
b = max(ax,cx);











