function out = logLnorm(obj)
     if isempty(obj.List), return; end
     lnorm = 0;
     for i=1:1:obj.nrTM
         lnorm = lnorm + logLnorm(obj.List{i});
     end
     if nargout == 1, out = lnorm;return;end
     obj.LnormEvaluation = lnorm;
end