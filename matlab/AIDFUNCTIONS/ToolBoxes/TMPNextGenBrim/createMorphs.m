function [FminMorphs,RegMorphs] = createMorphs(M,Pop,var,IndVar,BF,shape)
          FminMorphs = cell(2,var.nr2Boot);
          RegMorphs = cell(2,var.nr2Boot);
          RIP = updateRIP(Pop,M)';
          IndVar(:,var.Booting) = RIP;
          Mreg = getRegression(IndVar,Pop,(1:size(IndVar,2)));
          avgIndVar = mean(IndVar);
          avgY = mean(Pop);
          for i=1:1:var.nr2Boot
              avgRIP = mean(RIP(:,i));
              stdRIP = std(RIP(:,i));
              Mval = avgRIP-BF(i)*stdRIP;
              Pval = avgRIP+BF(i)*stdRIP;
              % fmin version
              C = getMorph(M(i,:),Mval);
              FminMorphs{1,i} = getScan(shape,C);
              C = getMorph(M(i,:),Pval);
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
%           for i=1:1:var.nr2Boot
%               for j=1:1:2
%                   viewer(FminMorphs{j,i});
%               end
%           end
%           
%           for i=1:1:var.nr2Boot
%               for j=1:1:2
%                   viewer(RegMorphs{j,i});
%               end
%           end
end

function [out,optX] = getMorph(M,val)      
         optX = fminsearch(@(X) MorphError(X,M,val),val);
         out = optX*M;
end

function out = MorphError(X,M,val)        
         C = X*M;
         test = updateRIP(C,M);
         out = abs(test-val);
end