function out = eval(obj,Tmodel)
         pdf(obj,Tmodel);
         [B,sumB] =  getB(obj,Tmodel); %#ok<NASGU>
         tmpout = -1*sum((1-B).*log(obj.PdfEvaluation));
         if nargout == 1, out = tmpout; return; end
         obj.Evaluation = tmpout;
end