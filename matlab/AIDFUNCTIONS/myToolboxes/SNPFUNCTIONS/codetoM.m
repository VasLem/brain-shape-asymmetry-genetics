function [out,flipped] = codetoM(geno)
         n0 = sum(geno==0);
         n2 = sum(geno==2);
         flipped = false;
         if n0<=n2
             out = geno;
         else
             out = abs(geno-2);
             out(out==3)=-1;
             flipped = true;
         end
end