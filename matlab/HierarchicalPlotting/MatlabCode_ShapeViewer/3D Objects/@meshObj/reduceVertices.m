function out = reduceVertices(obj,nrVertices,varargin)
    if nargout == 1
       obj = clone(obj);
       %obj.Visible = false;
       out = obj;
    end
    [Vindex,Findex] = getVindexFindex(obj,varargin{:}); %#ok<NASGU>
    if nrVertices >= length(Vindex), return; end
    factor = floor(obj.nrV/nrVertices);
    reduceTriangles(obj,floor(obj.nrF/factor));
end