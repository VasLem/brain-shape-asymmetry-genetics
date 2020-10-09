function [G1out,G2out] = ADASYN(G1,G2,DB,K)       
            [nr1,dim] = size(G1);
            [nr2,~] = size(G2);
            nrm = min(nr1,nr2);
            nrM = max(nr1,nr2);
            K = min(K,nrm-1);
            dcrit = nrm/nrM;
            if dcrit >= DB, G1out = G1;G2out = G2; return; end
            if nrM == nr1
               GM = G1;Gm = G2;
            else
               GM = G2;Gm = G1;
            end
            G = [Gm;GM];
            indM = (nrm+1:size(G,2));
            nrSyn = (nrM-nrm)*DB;
            distances = pdist2(Gm,G);
            [~,index] = sort(distances','ascend'); %#ok<UDIM>
            index = index(2:K+1,:);
            r = ismember(index,indM);
            r = sum(r)/K;
            if sum(r) > 0.01*nrSyn
                r = r/sum(r);
                g = round(r*nrSyn);
                totg = sum(g);
                Gadd = zeros(totg,dim);
                counter = 0;
                distances = squareform(pdist(Gm));
                [~,index] = sort(distances,'ascend');
                index = index(2:K+1,:);
                for i=1:1:nrm
                    if g(i)==0, continue; end
                    for l=1:g(i)
                        counter = counter+1;
                        nb = randi(K,1,1);
                        alpha = rand(1);
                        Gadd(counter,:) = Gm(i,:) + (Gm(index(nb,i),:)-Gm(i,:))*alpha;
                    end
                end
            else % close to seperate groups
                nrSyn = round(nrSyn);
                s = randsample(nrm,nrSyn,true);
                Gadd = zeros(nrSyn,dim);
                for f=1:1:nrSyn% generating new faces
                      nb = randsample(setdiff(1:nrm,s(f)),1);
                      alpha = rand(1);% random interpolation factor 
                      Gadd(f,:) = Gm(s(f),:) + (Gm(nb,:)-Gm(s(f),:))*alpha;% create new face
                end
            end
            Gm = [Gm; Gadd];
            if nrM == nr1
               G1out = GM;G2out = Gm;
            else
               G1out = Gm;G2out = GM;
            end
end