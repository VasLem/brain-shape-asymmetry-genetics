function out = createGeneEffects(A,shape,RefScan,BF,CondValues)
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
            [LocalE,LocalS] = getLocalStatistic(MR,stats.Yresiduals,E,'shape');
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
                eval(['morph' num2str(counter) ' = clone(scan);']);
            end
            out.Morph1 = morph1;
            out.Morph2 = morph2;
            out.STATS = STATS;
end

