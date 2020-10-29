function out = getProcrustesDirection(S1,S2)
         out = S2-S1;
         out = out/norm(out); 
end