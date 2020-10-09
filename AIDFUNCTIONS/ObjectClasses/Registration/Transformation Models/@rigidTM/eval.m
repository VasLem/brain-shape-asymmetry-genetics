function out = eval(obj,p)
    p = TM.convertInput(p);
    if isempty(obj.c),obj.c = p;end
    c = repmat(obj.c,1,p.nrV);    
    %out = rotation(obj)*(a-c)+c+translation(obj,a);
    p.Vertices = obj.Rotation*(p.Vertices-c)+c+repmat(obj.Translation,1,p.nrV);
    if nargout==1,out = p;return;end
    obj.Evaluation = p;% performs a copy;
    delete(p);% hence p needs to be deleted
end