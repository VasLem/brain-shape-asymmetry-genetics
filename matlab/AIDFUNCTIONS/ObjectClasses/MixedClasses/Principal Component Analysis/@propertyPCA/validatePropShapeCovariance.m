function out = validatePropShapeCovariance(obj,prop)

    Model = clone(obj.Model);
    out.trueprop = nan*zeros(1,obj.n);
    out.estprop = nan*zeros(1,obj.n);
    for i=1:1:obj.n
        tmp = Coeff2Struc(obj,obj.Tcoeff(1,:)');
        out.trueprop(i) = tmp.Prop(prop);
        clear tmp;
        
        
    end
    


end