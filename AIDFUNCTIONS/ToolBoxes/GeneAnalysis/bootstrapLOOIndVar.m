function [BIndVar,Iter] = bootstrapLOOIndVar(IndVar,DepVar,Model,type,index)
         BIndVar = IndVar;
         BIndVarNew = BIndVar;
         nrBoot = length(index);
         if nargout ==2, Iter(1).IndVar = BIndVarNew;end
         for bootstep=1:1:10
             disp(['Boot Step: ' num2str(bootstep)]);
             for i=1:1:nrBoot
                 BIndVarNew(:,index(i)) = myregressLOO(BIndVar,DepVar,Model,type,index(i))';
             end
             if bootstep == 1
                BIndVar = BIndVarNew;
                if nargout==2, Iter(1+bootstep).IndVar = BIndVarNew; end %#ok<*AGROW>
             else
                out = RV(BIndVar,BIndVarNew);
                disp(['RV: ' num2str(out.RV)]);
                BIndVar = BIndVarNew;
                if nargout==2, Iter(1+bootstep).IndVar = BIndVarNew; end
                if out.RV>=0.99
                   break;
                end
             end
         end
end