function out = updateCorrespondence(obj,Tmodel)
         % TM needs to be evaluated before SM eval !!!!
           %if isempty(Tmodel.Evaluation), error('TModel not evaluated, cannot perform correspondence update');end
           corresp = fastClone(Tmodel.Evaluation);
           corresp.Vertices = Tmodel.Evaluation.Vertices-repmat(obj.Evaluation,3,1).*obj.Direction;
           %corresp.Vertices = obj.TargetInfo(:,obj.Index);
           if nargout == 1, out = corresp; return; end
           obj.Correspondence = corresp;% performs a copy;
           delete(corresp);% hence corresp needs to be deleted
end