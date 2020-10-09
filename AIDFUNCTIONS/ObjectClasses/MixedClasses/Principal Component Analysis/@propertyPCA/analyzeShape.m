function out = analyzeShape(obj,Shape)
         out.ShapeCoeff = Struc2Coeff(obj.Model,Shape);
         vec = [out.ShapeCoeff; nan*ones(obj.nrP,1)];
         index = find(~isnan(vec));
         out.PropCoeff = weightedFit2(obj,vec,[],0.0001,index);
         out.Scan =getScan(obj,out.PropCoeff);
         maxp = getPlausibilty(obj.Model,obj.Model.AvgCoeff);
         p = getPlausibilty(obj.Model,out.ShapeCoeff);
         out.PercPlausible = (p/maxp)*100;
         out.Plausible = p;
end