function out = derivate(obj,Tmodel) 
    if isempty(obj.Smeasure), return; end
    if isempty(Tmodel.Evaluation),   error('TModel not evaluated');end
    for i=1:1:obj.nrSM
        derivate(obj.Smeasure{i},Tmodel);
    end
    if nargout == 1, out = obj.Derivative; return; end
end