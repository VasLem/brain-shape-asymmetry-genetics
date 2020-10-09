function [SP] = extractSymmetryPlane(obj)

    centerVertices(obj);
    cobj = clone(obj);
    cobj.Vertices(1,:) = -1*cobj.Vertices(1,:);
    T = rigidTM;
    match(T,obj,cobj);
    transform(T,cobj);


end