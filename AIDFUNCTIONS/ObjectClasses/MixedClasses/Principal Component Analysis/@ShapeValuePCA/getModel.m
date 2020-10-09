function out = getModel(obj,D)
    if nargout == 1,obj = clone(obj);out = obj;end
    if nargin==2&&isstruct(D)
       if isfield(D,'Shape')
            checkShape(obj);
            getModel(obj.Shape,D);
       end
       if isfield(D,'Value')
          checkValue(obj);
          getModel(obj.Value,D);
       end
    end
    D = obj.WS*eye(obj.nrSC)*(obj.Shape.Tcoeff');
    %D = obj.Shape.Tcoeff';
    D = [D;obj.Value.Tcoeff'];
    getModel@PCA(obj,D);
end