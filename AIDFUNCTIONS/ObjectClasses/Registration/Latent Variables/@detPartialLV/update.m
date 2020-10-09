function out = update(obj,Tmodel)
         if isempty(obj.CompleteP), error('Cannot update detOverLV:No CompleteP defined'); end
         if isempty(Tmodel.Evaluation), error('Cannot update detOverLV:Tmodel not evaluated'); end
         measures = obj.Function;    
         nrM = length(measures);
         tmpout = ones(1,nrM);
         nrReject = round(obj.TH*nrM);
         [Sm,Index] = sort(measures);
         Index = Index(end-nrReject:end);
         tmpout(Index) = 0;
         if nargout == 1, out = tmpout; return; end
         obj.Value = tmpout;
end