function out = estimatePropFromShapeCoeff(obj,ShapeCoeff)


         vec = [ShapeCoeff; nan*ones(obj.nrP,1)];
         index = find(~isnan(vec));
         
         
         out = obj.EigVec(index,:)'*(vec(index)-obj.AvgVec(index));
         %out = obj.EigVec(:,index)*(vec(index)-obj.AvgVec(index));
         
         out2 = obj.AvgVec + obj.EigVec*out;
         
         Scan = getScan(obj,out);
         
         propcoeff = Vec2Coeff(obj,vec,index);
         tmp = Coeff2Struc(obj,propcoeff);
         Scan = getScan(obj,propcoeff);
         out.Prop = Scan.UserData.Prop;
         out.propcoeff = propcoeff;
         %out = tmp.Prop;
%          out.propvec = Coeff2Vec(obj,propcoeff);
%          out.modelcoeff = vec;
%          Scan.UserData.Prop = out.Prop;
%          delete(tmp);


end