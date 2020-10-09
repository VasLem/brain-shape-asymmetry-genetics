function myWolfe(obj)
% Uses matlab opitmization toolbox linesearch stored in Utilities
[obj.S,f_new,fPrime_new,g_new,LSexitFlag,LSiter]=...
            matlabLineSearch({'fungrad',[],@ErrorFunction},...
            obj.X,obj.nrP,1,obj.nrP,obj.D,obj.F,obj.GtD,obj.S,obj.c1,obj.c2,-inf,100,...
            obj.ChangeTol,[],[],[],obj); %#ok<NASGU>
end

function [error,grad] = ErrorFunction(X,obj)
            obj.X = X;
            if nargout < 2
               fun(obj,'f');
               error = obj.F;
            else
               fun(obj,'fg');
               error = obj.F;
               grad = obj.G;
            end
end 