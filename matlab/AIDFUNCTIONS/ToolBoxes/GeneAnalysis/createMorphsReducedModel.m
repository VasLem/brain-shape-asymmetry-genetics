function [morphs,STATS] = createMorphsReducedModel(A,shape,RefScan,BF,CondValues)
         % multiple regression to obtain multiple R-Squared
            [~,~,~,~,~,pctvar] = plsregress(A,shape,size(A,2));
            R2Tot = sum(pctvar(2,:));
            STATS.R2Tot = R2Tot;
         % Conditioning regression to obtain reference face according to CondValues
            AC = A(:,1:end-1);
            [~,~,~,~,MC,pctvarC,~,statsC] = plsregress(AC,shape,size(AC,2));
            E = statsC.Yresiduals;
            mC = mean(AC);
            mShape = mean(shape);
            if isempty(CondValues)
               CondValues = mC; 
            end
            DX = mC-CondValues;
            DY = DX*MC(2:end,:);
            NY = mShape-DY;
            scanC = clone(RefScan);
            delete(scanC.PoseLM);
            scanC.Vertices = reshape(NY,3,length(NY)/3);
            R2Red = sum(pctvarC(2,:));
            R2Add = R2Tot-R2Red;
            STATS.R2Add = R2Add;
         % reduced model construction   
            AR = A(:,end);
            [~,~,~,~,~,~,~,statsR] = plsregress(AC,AR,1);
            AR = statsR.Yresiduals;
            [~,~,~,~,MR,pctvar,~,stats] = plsregress(AR,E,1);
            R2Part = pctvar(2);
            STATS.R2Part = R2Part;
            [LocalE,LocalS] = getLocalStatistic(MR,stats.Yresiduals,E);
            STATS.LocalE = LocalE;
            STATS.LocalS = LocalS;
         % creating morphs in reduced model space   
            Y = mean(E);
            m = mean(AR);mv = std(AR);
            range(1) = BF*(m-3*mv);
            range(2) = BF*(m+3*mv);
            counter = 0;
            for k=1:1:2
                V = range(k);
                counter = counter+1;
                DX = m-V;
                DY = DX*MR(2:end,:);
                NY = Y-DY;
                scan = clone(RefScan);
                delete(scan.PoseLM);
                scan.Vertices = scanC.Vertices + reshape(NY,3,length(NY)/3);
                morphs(counter).scan = clone(scan);
                if ~isempty(RefScan.PoseLM)
                    transferPoseLM(morphs(counter).scan,RefScan);
                end
            end
end

function [LocalE,LocalS] = getLocalStatistic(M,R,B)
         [n,nB] = size(B); 
         P = B-R;
         A = repmat(mean(B),n,1);
         SST = sum((B-A).^2);
         SSR = sum((P-A).^2);
         LocalE = zeros(4,nB/3);
         LocalE(1:3,:) = reshape(M(2,:),3,nB/3);
         LocalE(4,:) = sqrt(sum(LocalE(1:3,:).^2));
         SSR = reshape(SSR,3,nB/3);
         SSR = sum(SSR);
         SST = reshape(SST,3,nB/3);
         SST = sum(SST);
         LocalS = SSR./SST;
end