function [M,R2] = getBootSampleRegression(IndVar,DepVar,var,runs)
        R2out = true;
        nrS = size(IndVar,1);
        if nargout < 2, R2 = []; R2out = false;end
        for s=1:1:runs
            if runs==1% no random sampling
               SampleFi = 1:nrS; 
            else % random sampling with replacement
               SampleFi = randsample(nrS,nrS,true);
            end
            if s==1
               if R2out
                  [M,R2] = getRegression(IndVar(SampleFi,:),DepVar(SampleFi,:),var.Booting);
               else
                  M = getRegression(IndVar(SampleFi,:),DepVar(SampleFi,:),var.Booting);
               end
            else
               if R2out
                  [Mtmp,R2tmp] = getRegression(IndVar(SampleFi,:),DepVar(SampleFi,:),var.Booting);
                  R2 = R2+R2tmp;
               else
                  Mtmp = getRegression(IndVar(SampleFi,:),DepVar(SampleFi,:),var.Booting);
               end
               M = M+Mtmp;
            end
        end
        M = M/runs;
        R2 = R2/runs;
end