function out = pdf(obj,Tmodel)
         if isempty(obj.Smeasure), error('Cannot eval IP pdf, no Similarity Measure set'); end
         if ~(obj.Smeasure.nrV==Tmodel.nrV), eval(obj.Smeasure,Tmodel); disp('Smeasure evaluated');end
         %f = mvnpdf(obj.Smeasure.Evaluation',obj.Mu,obj.Sigma); Strange
         %mvnpdf of dim =1 is not equal to normpdf!!!!!
         f = ones(1,obj.nrV);
         for i=1:1:obj.nrSM
             f = f.*normpdf(obj.Smeasure.Evaluation(i,:),0,obj.PdfP(i));
         end
         if nargout == 1, out = f; return; end
         obj.PdfEvaluation = f;
end