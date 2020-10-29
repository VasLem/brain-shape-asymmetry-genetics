function out =  FallbackTest(obj)
    tmp = floor(obj.ObjFun.Evaluation*(1/obj.ChangeTol));
    out = ismember(tmp,obj.Old.Evals);
    obj.Old.Evals(obj.MstepEval + 1) = tmp;
end