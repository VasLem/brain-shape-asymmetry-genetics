function [FminMorphs,RegMorphs] = extractMorphs(M,Pop,var,IndVar,BF,shape)
%[FminMorphs,RegMorphs] = extractMorphs(BRIM.HMedM,BRIM.DepVar,BRIM.Var,BRIM.IndVar,1,SymSpace);
          FminMorphs = cell(2,var.nr2Boot);
          RegMorphs = cell(2,var.nr2Boot);
          RIP = updateRIP(Pop,M)';
          %RIP = CRIP;
          disp('OK');
          IndVar(:,var.Booting) = RIP;
          %PopN = Pop.*repmat(shape.EigStd',size(Pop,1),1);
          Mreg = getRegression(IndVar,Pop,(1:size(IndVar,2)));
          avgIndVar = nanmean(IndVar);
          avgY = nanmean(Pop);
          for i=1:1:var.nr2Boot
              avgRIP = nanmean(RIP(:,i));
              stdRIP = nanstd(RIP(:,i));
              Mval = avgRIP-BF(i)*stdRIP;
              Pval = avgRIP+BF(i)*stdRIP;
              % fmin version
              C = getFminMorph(M(i,:),Mval);
              FminMorphs{1,i} = getScan(shape,C);
              C = getFminMorph(M(i,:),Pval);
              FminMorphs{2,i} = getScan(shape,C);
              % regression version
              X = avgIndVar;
              X(var.Booting(i)) = Mval;
              deltaX = X-avgIndVar;
              dY = deltaX*Mreg;
              C = avgY+dY;
              RegMorphs{1,i} = getScan(shape,C);
              X(var.Booting(i)) = Pval;
              deltaX = X-avgIndVar;
              dY = deltaX*Mreg;
              C = avgY+dY;
              RegMorphs{2,i} = getScan(shape,C);
          end
end

function [out,optX] = getFminMorph(M,val)      
         optX = fminsearch(@(X) MorphError(X,M,val),val);
         out = optX*M;
end

function out = MorphError(X,M,val)        
         C = X*M;
         test = updateRIP(C,M);
         out = abs(test-val);
end