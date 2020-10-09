function out = eval(obj,Tmodel) 
    if obj.nrSM == 0, return; end
    if isempty(Tmodel.Evaluation),   error('TModel not evaluated');end
    for i=1:1:obj.nrSM
        eval(obj.Smeasure{i},Tmodel);
    end
    if nargout == 1, out = obj.Evaluation; return; end
end