function out = update(obj,Tmodel)
         if isempty(obj.CompleteP), error('Cannot update detProjectionLV:No CompleteP defined'); end
         if isempty(Tmodel.Evaluation), error('Cannot update detProjectionLV:Tmodel not evaluated'); end
         tmpout = ones(1,Tmodel.nrV);
         % first projection test
            angle = angleGradientViewDirection(obj.Function,Tmodel.Evaluation);
            tmpout((angle>=85)) = 0;
         % second projection test
            val = getValues(obj.Function,Tmodel.Evaluation);
            tmpout(isnan(val(1,:))) = 0;
         if nargout == 1, out = tmpout; return; end
         obj.Value = tmpout;
end