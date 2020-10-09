function out = update(obj,Tmodel)
         if isempty(obj.CompleteP), error('Cannot update detOverLV:No CompleteP defined'); end
         if isempty(Tmodel.Evaluation), error('Cannot update detOverLV:Tmodel not evaluated'); end
         Index = KNN(obj.Function,Tmodel.Evaluation,1);    
         %Index = knn(kde(obj.Function.Vertices,5),Tmodel.Evaluation.Vertices,1);
         tmpout = double(~ismember(Index,obj.Function.Border.VerticesIndex));
         if nargout == 1, out = tmpout; return; end
         obj.Value = tmpout;
end