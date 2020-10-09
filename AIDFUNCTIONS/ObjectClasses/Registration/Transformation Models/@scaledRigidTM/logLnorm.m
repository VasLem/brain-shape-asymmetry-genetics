function out = logLnorm(obj)
     lnorm = 0;
     if nargout == 1, out = lnorm;return;end
     obj.LnormEvaluation = lnorm;
end