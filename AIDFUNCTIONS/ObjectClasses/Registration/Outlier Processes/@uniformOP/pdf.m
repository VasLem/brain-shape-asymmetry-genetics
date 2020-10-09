function out = pdf(obj,Tmodel)
         f = repmat(obj.Level,1,Tmodel.nrV);
         if nargout == 1, out = f; return; end
         obj.PdfEvaluation = f;
end