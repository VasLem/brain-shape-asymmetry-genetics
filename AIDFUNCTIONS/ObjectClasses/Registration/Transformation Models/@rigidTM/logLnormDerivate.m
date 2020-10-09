function out = logLnormDerivate(obj)
     lnormgrad = 0;
     if nargout == 1, out = lnormgrad;return; end
     obj.LnormDerivative = lnormgrad;
     % function obj.ActiveP
end