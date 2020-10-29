function [G1out,G2out] = DownSample(G1,G2,DB)       
            nr1 = size(G1,1);
            nr2 = size(G2,1);
            nrm = min(nr1,nr2);
            nrM = max(nr1,nr2);
            dcrit = nrm/nrM;
            if dcrit >= DB, G1out = G1;G2out = G2; return; end
            nrKeep = round(nrm/DB);
            s = randsample(nrM,nrKeep,false);
            if nrM == nr1
               G1out = G1(s,:);G2out = G2;
            else
               G1out = G1;G2out = G2(s,:);
            end
end