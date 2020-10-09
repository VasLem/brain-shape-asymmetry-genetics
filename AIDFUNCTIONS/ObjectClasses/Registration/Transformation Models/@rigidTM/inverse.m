function out = inverse(obj,p)
    p = TM.convertInput(p);
    %if isempty(obj.c),obj.c = p;end
    c = repmat(obj.c,1,p.nrV);
    Rotation = obj.Rotation^-1;
    Translation = obj.Translation;
    p.Vertices = Rotation*(p.Vertices-c-repmat(Translation,1,p.nrV))+c;
    out = p;
end