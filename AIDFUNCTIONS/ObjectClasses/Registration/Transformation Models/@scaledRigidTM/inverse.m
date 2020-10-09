function out = inverse(obj,p)
    p = TM.convertInput(p);
    invRotation = obj.Rotation^-1;
    p.Vertices = invRotation*(p.Vertices/obj.Scale-repmat(obj.Translation,1,p.nrV));
    out = p;
end