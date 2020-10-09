function [obj,NIndVar,Iter] = BRIMC(obj)
         if ~checkNrCObservations(obj), error('different amount of observations between C and Y'); end
         [NIndVar,Iter] = BRIM(obj,obj.C,obj.Y,(1:obj.nrC),[]);
         obj.RIPC = NIndVar;
end