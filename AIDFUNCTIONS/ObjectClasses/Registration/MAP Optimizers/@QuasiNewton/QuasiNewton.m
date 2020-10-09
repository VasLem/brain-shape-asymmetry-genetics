classdef QuasiNewton < Gradient
    % This is the abstract interface class for MAP function optimizers
    % based on a Gradient Steepest Descent during Mstep of
    % the optimization
   properties
       HessianInitial = 'scaled';% What is the initial estimate 'hessian', 'eye', 'scaled'
       HessianMethod = 'BFGS';% Method to update Hessian estimation
       HessianCorr = 100;% Number of Hessian corrections to store in memory (relevant for LBFGS)
       Hdiag = [];% needed for LBFGS
       R = [];% Current Hessian estimation
       Damped = true;% use damped update , only relevant for BFGS
   end
   properties (Dependent = true)
       BestStep;
       BoundStep;
   end
   methods %Constructor
        function obj = QuasiNewton(varargin)
          obj = obj@Gradient(varargin{:});
        end
   end
   methods % Special Getting & Setting
       function out = get.BestStep(obj) %#ok<INUSD>
           switch obj.HessianMethod
               case 'BFGS'
                   switch obj.HessianInitial
                       case 'hessian'
                           out = 1;
                       otherwise
                           out = min(1,1/sum(abs(obj.G)));
                   end
               case 'LBFGS'
                   out = min(1,1/sum(abs(obj.G)));
           end
       end
       function out = get.BoundStep(obj)
          %out = obj.S;% No bounding required
          out = min(obj.S,1e4/(1+sum(abs(obj.G))));
       end
       function obj = set.HessianMethod(obj,in)
           if ~strcmp(in,'BFGS')&&...
              ~strcmp(in,'LBFGS')
               error('Wrong Type of HessianMethod, BFGS or LBFGS is required');
           else
              obj.HessianMethod = in;
           end
       end
       function obj = set.HessianInitial(obj,in)
           if ~strcmp(in,'hessian')&&...
              ~strcmp(in,'eye')&&...
              ~strcmp(in,'scaled')
               error('Wrong Type of HessianInitial, hessian or eye is required');
           else
              obj.HessianInitial = in;
           end
       end
   end
   methods % InterFace functions
       function descentDirection(obj)
           switch obj.HessianMethod
               case 'BFGS'
                   bfgsDirection(obj);
               case 'LBFGS'
                   lbfgsDirection(obj);
           end
       end
       function bfgsDirection(obj)
           if obj.MstepEval == 0
              % set initial hessian for hessian
              switch obj.HessianInitial
                  case 'hessian'
                     % Update Hessian
                       fun(obj, 'h');
                     % Modified Incomplete Cholesky
                       obj.R = mcholinc(obj.H);
                       obj.D = -obj.R\(obj.R'\obj.G);
                  otherwise
                     % use steepest descent
                       obj.D = -obj.G;
              end
              obj.Old.G = obj.G;
              return;
           end
           y = obj.G-obj.Old.G;
           if isempty(obj.S), obj.S = 0.001; end
           s = obj.S*obj.D;
           if obj.MstepEval == 1
               % set initial Hessian approximation for eye and scaled
               switch obj.HessianInitial
                  case 'eye'
                      obj.R = eye(length(obj.G));
                   case 'scaled'
                      obj.R = sqrt((y'*y)/(y'*s))*eye(length(obj.G));
                   otherwise
                      % do nothing (hessian already done)
               end             
           end
           if obj.Damped
               eta = .02;
               B = obj.R'*obj.R;
               if y'*s < eta*s'*B*s
                 theta = min(max(0,((1-eta)*s'*B*s)/(s'*B*s - y'*s)),1);
                 y = theta*y + (1-theta)*B*s;
               end
               [obj.R,posDef]= cholupdate(cholupdate(obj.R,y/sqrt(y'*s)),obj.R'*obj.R*s/sqrt(s'*obj.R'*obj.R*s),'-');%#ok<NASGU>
           else
               if y'*s > 1e-10
                  [obj.R,posDef]= cholupdate(cholupdate(obj.R,y/sqrt(y'*s)),obj.R'*obj.R*s/sqrt(s'*obj.R'*obj.R*s),'-');%#ok<NASGU>
               else
                   % skipping update
               end
           end
           %if posDef
            obj.D = -obj.R\(obj.R'\obj.G);
%            else
%             obj.D = -obj.G;
%            end
           obj.Old.G = obj.G;
       end
       function lbfgsDirection(obj)
           if obj.MstepEval < 20
              obj.ChangeTol = 0.01;
           else
              obj.ChangeTol = 1;
           end
           if obj.MstepEval == 0
              obj.D = -obj.G;
              obj.Old.dirs = zeros(length(obj.G),0);
              obj.Old.stps = zeros(length(obj.D),0);
              obj.Hdiag = 1;
           else
              [obj.Old.dirs,obj.Old.stps,obj.Hdiag] = lbfgsUpdate(obj.G-obj.Old.G,obj.S*obj.D,obj.HessianCorr,false,obj.Old.dirs,obj.Old.stps,obj.Hdiag);
              obj.D = lbfgs(-obj.G,obj.Old.dirs,obj.Old.stps,obj.Hdiag);
              %obj.D = lbfgsC(-obj.G,obj.Old.dirs,obj.Old.stps,obj.Hdiag);
           end
           obj.Old.G = obj.G;
              
       end
       function defaultOptions(obj)
           obj.LS = 'bracketing';
           obj.Bracketing = 'wolfe';
           obj.LSinit = 'newton';
           obj.HessianInitial = 'scaled';
           obj.Damped = true;
           %obj.ChangeTol = 1;
           if obj.nrP < 1000
              obj.HessianMethod = 'BFGS';    
           else
              obj.HessianMethod = 'LBFGS';
           end
           obj.c2 = 0.9;
           %obj.MaxIter = 100;    
       end 
       function copy(obj,cobj)
           copy@Gradient(obj,cobj);
           cobj.HessianInitial = obj.HessianInitial;
           cobj.HessianMethod = obj.HessianMethod;
           cobj.HessianCorr = obj.HessianCorr;
           cobj.Hdiag = obj.Hdiag;
           cobj.R = obj.R;
           cobj.Damped = obj.Damped;
       end
   end   
end % classdef