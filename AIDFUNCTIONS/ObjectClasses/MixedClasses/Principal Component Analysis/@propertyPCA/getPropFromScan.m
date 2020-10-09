function [out,propcoeff] = getPropFromScan(obj,Scan)        
         T = rigidTM;
         match(T,obj.Model.Average,Scan);
         tmp = transform(T,Scan);
         vec = getCoeff(obj.Model,tmp);
         vec = [vec; nan*ones(obj.nrP,1)];
         propcoeff = partialFit(obj,vec);
         out = Coeff2Struc(obj,propcoeff);
         out.propvec = Coeff2Vec(obj,propcoeff);
         out.modelcoeff = vec;
         Scan.UserData.Prop = out.Prop;
         delete(tmp);
end