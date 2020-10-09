function out = update(obj,Tmodel)
         if isempty(obj.CompleteP)||isempty(obj.CompleteP.InlierP)||isempty(obj.CompleteP.OutlierP), return; end
         pdf(obj.CompleteP.InlierP,Tmodel);
         pdf(obj.CompleteP.OutlierP,Tmodel);        
         tmpout = 1-exp(-1*exp(-1*(obj.CompleteP.InlierP.Smeasure.Evaluation/obj.CompleteP.InlierP.Sigma)));
         if nargout == 1, out = tmpout; return; end
         obj.Value = tmpout;
end