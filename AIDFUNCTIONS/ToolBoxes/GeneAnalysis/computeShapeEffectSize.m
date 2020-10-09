function out = computeShapeEffectSize(A,shape)
         % multiple regression to obtain multiple R-Squared
            [~,~,~,~,~,pctvar] = plsregress(A,shape,size(A,2));
            R2Tot = sum(pctvar(2,:));
         % reduced model computation (it is always the last one to test)
            AR = A(:,1:end-1);
            [~,~,~,~,~,pctvar,~,statsR] = plsregress(AR,shape,size(AR,2));
            E = statsR.Yresiduals;
            R2Red = sum(pctvar(2,:));
            R2Add = R2Tot-R2Red;


end