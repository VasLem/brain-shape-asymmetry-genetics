function out = update(obj,Tmodel)
         if isempty(obj.CompleteP)||isempty(obj.CompleteP.InlierP)||isempty(obj.CompleteP.OutlierP), return; end
         pdf(obj.CompleteP.InlierP,Tmodel);
         pdf(obj.CompleteP.OutlierP,Tmodel);
         fpdf = obj.CompleteP.InlierP.PdfEvaluation;
         gpdf = obj.CompleteP.OutlierP.PdfEvaluation;
         tmpout = fpdf./(fpdf+gpdf);
         if nargout == 1, out = tmpout; return; end
         obj.Value = tmpout;
end