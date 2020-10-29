function ArmijoBacktrack(obj)
%
% Backtracking linesearch to satisfy Armijo condition
%
% Inputs:
%   X: starting location
%   S: initial step size
%   D: descent direction
%   F: function value at starting location
%   Fref: reference function value (usually funObj(x))
%   GtD: directional derivative at starting location
%   c1: sufficient decrease parameter
%   debug: display debugging information
%   Backtrack: type of interpolation
%   ChangeTol: minimum allowable step length
%   Plot: do a graphical display of interpolation
%   ObjFun: objective function
%
% Outputs:
%   S: step length
%   F_new: function value at x+t*d
%   G_new: gradient value at x+t*d
%   H: Hessian at initial guess (only computed if requested

    if any(isnan(obj.D))
        return;
    end
    f = obj.F;
    g = obj.G;
    dr = obj.D;
    x = obj.X;
    % Evaluate the Objective and Gradient at the Initial Step
    obj.X = x+obj.S*obj.D;
    if strcmp(obj.Backtrack, 'cubic fg')
        fun(obj,'fg');
    else
        fun(obj,'f');
    end
    funEvals = 1;
    while obj.F > obj.Fref + obj.c1*obj.S*obj.GtD || ~isLegal(obj.F)

        temp = obj.S;
        if strcmp(obj.Backtrack,'step halving') || ~isLegal(obj.F)
            % Backtrack w/ fixed backtracking rate
            obj.S = 0.5*obj.S;
        elseif strcmp(obj.Backtrack, 'cubic fg') && isLegal(obj.G)
            % Backtracking w/ cubic interpolation w/ derivative
            obj.S = polyinterp([0 f obj.GtD; obj.S obj.F obj.G'*dr],obj.Plot);
        elseif strcmp(obj.Backtrack,'cubic f')
          if funEvals < 2 || ~isLegal(f_prev)
                % Backtracking w/ quadratic interpolation (no derivative at
                % new point)
                obj.S = polyinterp([0 f obj.GtD; obj.S obj.F sqrt(-1)],obj.Plot);
          else
                % Backtracking w/ cubic interpolation (no derivatives at
                % new points)
                obj.S = polyinterp([0 f obj.GtD; obj.S obj.F sqrt(-1); s_prev f_prev sqrt(-1)],obj.Plot);
          end
        else
           error('wrong type of Backtracking');
        end

        % Adjust if change in t is too small/large

        if obj.S < temp*1e-3
            disp('Step to small');
            obj.S = temp*1e-3;
        elseif obj.S > temp*0.6
            disp('step to big');
            obj.S = temp*0.6;
        end

        f_prev = obj.F;
        s_prev = temp;
        obj.X = x+obj.S*obj.D;
        if strcmp(obj.Backtrack, 'cubic fg')
            fun(obj,'fg');
        else
            fun(obj,'f');
        end
        funEvals = funEvals+1;

        % Check whether step size has become too small
        if sum(abs(obj.S*obj.D)) <= obj.ChangeTol;
            disp('Backtracking Line Search Failed\n');
            obj.S = 0;
            obj.X = x;
            fun(obj,'f');
            break;
        end
    end
end

function [legal] = isLegal(v)
legal = sum(any(imag(v(:))))==0 & sum(isnan(v(:)))==0 & sum(isinf(v(:)))==0;
end