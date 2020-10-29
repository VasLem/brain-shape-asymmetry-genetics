function [out,optX] = getMorph(M,val)      
         optX = fminsearch(@(X) MorphError(X,M,val),val);
         out = optX*M;
end

function out = MorphError(X,M,val)        
         C = X*M;
         test = updateRIP(C,M);
         out = abs(test-val);
end