function out = update(obj,Tmodel) %#ok<INUSD>
         if isempty(obj.CompleteP) || isempty(obj.CompleteP.InlierP)
            L = obj.Level;% No update, possible
         else
            inlierP = obj.CompleteP.InlierP;
            switch inlierP.Type
                case 'gaussianIP'
                    L = (1/sqrt(((2*pi)^2)*det(inlierP.Sigma)))*exp(-0.5*inlierP.Kappa^2);
                case 'KdeIP'
                    L = ((1/sqrt(((2*pi)^2)*det(inlierP.Sigma)))*exp(-0.5*inlierP.Kappa^2))/inlierP.Smeasure.TargetInfo.Npts;
                otherwise
                    L = obj.Level;% Currently no update
            end
         end
         if nargout == 1, out = L; return; end
         obj.Level = L;
end