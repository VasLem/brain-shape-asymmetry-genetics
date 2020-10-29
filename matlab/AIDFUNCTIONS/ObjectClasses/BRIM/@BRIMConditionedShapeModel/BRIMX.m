function [obj,NIndVar,Iter] = BRIMX(obj,index)
         if ~checkNrXObservations(obj), error('different amount of observations between X and Y'); end    
         [NIndVar,Iter] = BRIM(obj,[obj.RIPC obj.X(:,index)],obj.Y,index,(1:obj.nrC));
         obj.RIPX(:,index) = NIndVar(:,obj.nrC+1:end);
end