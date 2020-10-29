function getAverage(obj,D)
         if ~isstruct(D), return; end
         if isfield(D,'Shape')
            checkShape(obj);
            getAverage(obj.Shape,D);
         end
         if isfield(D,'Value')
            checkValue(obj);
            getAverage(obj.Value,D);
         end
end