function out = eval(obj,Tmodel)
         % TM needs to be evaluated before SM eval !!!!
         % Target of SM needs to be set before SM eval!!!
         % if not [] values are returned
         if isempty(obj.Target), error('Target not set for SM'); end
         if isempty(Tmodel.Evaluation),   error('Cannot evaluate Smeasure, because TModel not evaluated');end         
         updateClosestPoints(obj,Tmodel);
         val = obj.Distance-Tmodel.Evaluation.Distance;
         if nargout == 1, out = val; end
         obj.Evaluation = val;       
end
