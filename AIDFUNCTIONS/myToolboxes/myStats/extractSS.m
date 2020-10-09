function out = extractSS(D,a,b,n)
         N = a*b*n;
         T = 0;A = 0;B = 0;R = 0;
         for i=1:1:(N-1)
             for j=(i+1):1:N
                 d = D(i,j)^2;
                 T = T+d;
                 if inA(i,N,a)==inA(j,N,a)
                    A = A+d;
                 end
                 if inB(i,N,a,b)==inB(j,N,a,b)
                    B = B+d;
                    if inA(i,N,a)==inA(j,N,a)
                       R = R+d;
                    end
                 end
             end
         end
         out.SST = T/N;
         out.SSAW = A/(n*b);
         out.SSA = out.SST-out.SSAW;
         out.SSBW = B/(n*a);
         out.SSB = out.SST-out.SSBW;
         out.SSR = R/n;
         out.SSAB = out.SST-out.SSA-out.SSB-out.SSR;
         out.MSA = out.SSA/(a-1);
         out.MSB = out.SSB/(b-1);
         out.MSAB = out.SSAB/((a-1)*(b-1));
         out.MSR = out.SSR/(N-a*b);
end
function out = inA(i,N,a)
         nA = N/a;
         out = ceil(i/nA);
end
function out = inB(i,N,a,b)
         nA = N/a;
         nB = (N/a)/b;
         i = i-((inA(i,N,a)-1)*nA);
         out = ceil(i/nB);
end
