classdef Gradient < mapOPT
    % This is the abstract interface class for MAP function optimizers
    % based on a Gradient Descent during Mstep of
    % the optimization
   properties
       D = [];% Direction, Dependent on method (see subclasses)
       S = [];% Actual Step
       GtD = [];% G transpose times D
       LS = 'none';% 'none' , 'backtrack', 'bracketing'
       Backtrack = 'step halving';% 'step halving', 'cubic f', 'cubic fg';
       Bracketing = 'wolfe'% 'fminbnd', 'wolfe'
       Plot = false; % do a graphical display of interpolation, backtrack or bracketing
       LSinit = 'newton'; % Initial stepsize for LS, 'newton', 'similar' (to previous), 'quadratic' (based on function value and new value gradient)
       Computation = 'Analytical'; % How to compute gradient, Analytical or Numeric;
       Fref = [];% Using non-monotone Armijo backtrack
       c1 = 1e-4;% Sufficient Decrease for Armijo (backtrack) condition (1e-4)
       c2 = 0.9;% Curvature Decrease for Wolfe conditions (.2 for conjugate Gradient, .9 otherwise)
   end
   properties (GetAccess = private, SetAccess = private)
       BestStep;% used to initialize step guess
       BoundStep;% used to bound step
   end
   properties (Dependent = true)
       nrP;% number of parameters
       F;% function evaluation
       G;% function Gradient
       H;% function Hessian
       X;% Current Parameters
   end
   methods %Constructor
        function obj = Gradient(varargin)
          obj = obj@mapOPT(varargin{:});
          if nargin>1
             Input = find(strcmp(varargin, 'D'));if ~isempty(Input), obj.D = varargin{Input+1}; end;
             Input = find(strcmp(varargin, 'S'));if ~isempty(Input), obj.S = varargin{Input+1}; end;
             Input = find(strcmp(varargin, 'Old'));if ~isempty(Input), obj.Old = varargin{Input+1}; end;
             Input = find(strcmp(varargin, 'Backtrack'));if ~isempty(Input), obj.Backtrack = varargin{Input+1}; end;
             Input = find(strcmp(varargin, 'Bracketing'));if ~isempty(Input), obj.Bracketing = varargin{Input+1}; end;
             Input = find(strcmp(varargin, 'Plot'));if ~isempty(Input), obj.Plot = varargin{Input+1}; end;
             Input = find(strcmp(varargin, 'LS'));if ~isempty(Input), obj.LS = varargin{Input+1}; end;
             Input = find(strcmp(varargin, 'LSinit'));if ~isempty(Input), obj.LSinit = varargin{Input+1}; end;
             Input = find(strcmp(varargin, 'Computation'));if ~isempty(Input), obj.Computation = varargin{Input+1}; end;
             
          end
        end
   end
   methods % Special Setting and Getting
       function out = get.nrP(obj)
           if isempty(obj.ObjFun), out = 0; return; end
           out = obj.ObjFun.nrP;
       end
       function out = get.F(obj)
           if isempty(obj.ObjFun), out = []; return; end
           out = obj.ObjFun.Evaluation;
       end
       function obj = set.F(obj,in)
           if isempty(obj.ObjFun), return; end
           obj.ObjFun.Evaluation = in;
       end
       function out = get.G(obj)
           if isempty(obj.ObjFun), out = []; return; end
           out = obj.ObjFun.Derivative;
       end
       function obj = set.G(obj,in)
           if isempty(obj.ObjFun), return; end
           obj.ObjFun.Derivative = in;
       end
       function out = get.H(obj)
           if isempty(obj.ObjFun), out = []; return; end
           out = obj.ObjFun.Hessian;
       end
       function obj = set.H(obj,in)
           if isempty(obj.ObjFun), return; end
           obj.ObjFun.Hessian = in;
       end
       function out = get.X(obj)
           if isempty(obj.ObjFun), out = []; return; end
           out = obj.ObjFun.P;
       end
       function obj = set.X(obj,in)
           if isempty(obj.ObjFun), return; end
           obj.ObjFun.P = in;
       end
       function obj = set.LS(obj,in)
           if ~strcmp(in,'none')&&...
              ~strcmp(in,'backtrack')&&...
              ~strcmp(in,'bracketing')
              error('Wrong Type of Line Search, none, backtrack or bracketing is required');
           else
              obj.LS = in;
           end
       end
       function obj = set.Backtrack(obj,in)
           if ~strcmp(in,'step halving')&&...
              ~strcmp(in,'cubic fg')&&...
              ~strcmp(in,'cubic f')
               error('Wrong Type of Backtrack, step halving, cubic fg or cubic f is required');
           else
              obj.Backtrack = in;
           end
       end
       function obj = set.Bracketing(obj,in)
           if ~strcmp(in,'wolfe')&&...
              ~strcmp(in,'fminbnd')
               error('Wrong Type of Bracketing, wolfe or fminbnd is required');
           else
              obj.Bracketing = in;
           end
       end
   end
   methods % InterFace functions
       function out = initialize(obj)
           if nargout == 1
                    obj = clone(obj);
                    out = obj;
           end
           initialize@mapOPT(obj);
           obj.S = [];
           obj.GtD = [];
           obj.D = zeros(obj.nrP,1);
           obj.Fref = [];
       end      
       function Mstep(obj)
           % ******* INITIALIZE **********
                initializeMStep(obj);
             % Check Optimality Condition
                if strcmp(obj.ChangeType, 'Evaluation')
                    if sum(abs(obj.G)) <= obj.ChangeTol
                       %disp('Optimality condition reached');
                       %obj.ExitFlag = [obj.ExitFlag 'Optimality condition reached \r'];
                       return;
                    end
                end
           % ******* GET DESCENT DIRECTION **********
                descentDirection(obj);% Subclass Routine
           % ******* COMPUTE STEP LENGTH **********
                obj.Old.F = obj.F;
                obj.GtD = obj.G'*obj.D;%Directional Derivative
             % Check that progress can be made along direction
                if strcmp(obj.ChangeType, 'Evaluation')
                    if obj.GtD > -obj.ChangeTol
                        %obj.ExitFlag = [obj.ExitFlag 'Directional Derivative below ChangeTol \r'];
                        return;
                    end
                end
              % Select initial Guess
                stepGuess(obj);% Gradient Class routine
              % perform Line Search
                lineSearch(obj);% Gradient Class routine
       end
       function initializeMStep(obj)
           %copy(MapFunction,obj.ObjFun);
           if isempty(obj.Fref)
              obj.Fref = obj.F;
           else
              obj.Fref = max(obj.Fref,obj.F);
           end
           fun(obj,'g');% Get function Gradient
       end
       function fun(obj,fgh)
              if ~isempty(strfind(fgh,'f'))
                 Ffun(obj);
              end
              if ~isempty(strfind(fgh,'g'))
                 Gfun(obj);
              end
              if ~isempty(strfind(fgh,'h'))
                 Hfun(obj);
              end
       end
       function Ffun(obj)
                % Function evaluation at X
                  eval(obj.ObjFun);
       end
       function Gfun(obj)
                % Function derivative at X
                switch obj.Computation
                   case 'Analytical'
                      derivate(obj.ObjFun);
                   case 'Numeric'
                      autoGrad(obj);
                   otherwise
                      error('wrong type of derivative computation')
                end
       end
       function Hfun(obj)
                % Function Hessian at X
                switch obj.Computation
                   case 'Analytical'
                      semiAutoHess(obj);
                   case 'Numeric'
                      autoHess(obj);
                   otherwise
                      error('wrong type of derivative computation')
                end
       end
       function autoGrad(obj)
           % requires Ffun first to obtain gradient in X
           f = obj.ObjFun.Evaluation;
           mu = 2*sqrt(1e-12)*(1+norm(obj.ObjFun.P))/norm(obj.ObjFun.nrP);
           diff = zeros(obj.ObjFun.nrP,1);
           M = clone(obj.ObjFun);
           for j = 1:1:obj.ObjFun.nrP
             ej = zeros(obj.ObjFun.nrP,1);
             ej(j) = 1;
             M.P = obj.ObjFun.P+mu*ej;
             diff(j,1) = eval(M);
           end
           obj.ObjFun.Derivative = (diff-f)/mu;
           delete(M);    
       end
       function autoHess(obj)
           % requires Ffun and Fgrad first
           g = obj.ObjFun.Derivative;
           mu = 2*sqrt(1e-12)*(1+norm(obj.ObjFun.P))/norm(obj.ObjFun.nrP);
           diff = zeros(obj.ObjFun.nrP);
           M = clone(obj.ObjFun);
           for j = 1:1:obj.ObjFun.nrP
               ej = zeros(obj.ObjFun.nrP,1);
               ej(j) = 1;
               M.P = obj.ObjFun.P+mu*ej;
               eval(M);
               Gradient.autoGrad(M);
               diff(:,j) = M.Derivative;
           end
           delete(M);
           H = (diff-repmat(g,[1 obj.ObjFun.nrP]))/mu;
           obj.ObjFun.Hessian = (H+H')/2;
       end
       function semiAutoHess(obj)
           % requires Ffun and Fgrad first
           g = obj.ObjFun.Derivative;
           mu = 2*sqrt(1e-12)*(1+norm(obj.ObjFun.P))/norm(obj.ObjFun.nrP);
           diff = zeros(obj.ObjFun.nrP);
           M = clone(obj.ObjFun);
           for j = 1:1:obj.ObjFun.nrP
               ej = zeros(obj.ObjFun.nrP,1);
               ej(j) = 1;
               M.P = obj.ObjFun.P+mu*ej;
               eval(M);
               derivate(M);
               diff(:,j) = M.Derivative;
           end
           delete(M);
           H = (diff-repmat(g,[1 obj.ObjFun.nrP]))/mu;
           obj.ObjFun.Hessian = (H+H')/2;
       end
       function str = strInfo(obj)
            str = strInfo@mapOPT(obj);
            str = [str '/Step: ' num2str(obj.S)];
       end
       function copy(obj,cobj)
           copy@mapOPT(obj,cobj);
           cobj.D = obj.D;
           cobj.S = obj.S;
           cobj.GtD = obj.GtD;
           cobj.LS = obj.LS;
           cobj.Backtrack = obj.Backtrack;
           cobj.Bracketing = obj.Bracketing;
           cobj.Plot = obj.Plot;
           cobj.LSinit = obj.LSinit;
           cobj.Computation = obj.Computation;
           cobj.Fref = obj.Fref;
           cobj.c1 = obj.c1;
           cobj.c2 = obj.c2;
       end
   end
   methods (Abstract = true)% Abstract Interface functions
       defaultOptions(obj);% Setting Default options         
       descentDirection(obj);% updating direction, dependent on method of concrete subclass
   end
end % classdef