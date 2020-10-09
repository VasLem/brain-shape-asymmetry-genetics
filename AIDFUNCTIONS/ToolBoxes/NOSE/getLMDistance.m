function out = getLMDistance(LM1,LM2)
         out = sqrt(sum((LM1-LM2).^2));
end