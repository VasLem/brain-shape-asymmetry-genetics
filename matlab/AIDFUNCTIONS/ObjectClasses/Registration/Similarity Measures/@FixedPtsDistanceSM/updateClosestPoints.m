function updateClosestPoints(obj,Tmodel)
    if ~isempty(obj.Index)
        obj.Difference = Tmodel.Evaluation.Vertices-obj.TargetInfo(:,obj.Index);
    else
        obj.Difference = Tmodel.Evaluation.Vertices-obj.TargetInfo;
    end
    obj.Distance = sqrt(sum(obj.Difference.^2));    
end