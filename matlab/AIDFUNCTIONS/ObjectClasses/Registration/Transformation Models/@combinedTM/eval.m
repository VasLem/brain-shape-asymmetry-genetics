function out = eval(obj,p)
    if obj.nrTM==0, return; end
    for i=obj.nrTM:-1:1% from back to front
        eval(obj.List{i},p);
        p = obj.List{i}.Evaluation;
    end
    if nargout==1,out=p;return;end
    obj.Evaluation = p;% performs a copy;
    %delete(p);% p does not need to be deleted, was not copyied in this
    %routine, only in subroutines!
end