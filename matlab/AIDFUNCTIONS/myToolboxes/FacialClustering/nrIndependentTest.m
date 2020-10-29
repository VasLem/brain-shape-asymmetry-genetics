function out = nrIndependentTest(COR)
         E = real(eig(COR));
         E = E(find(E>0));
         out  = sum((E>=1) + (E-floor(E)));
end