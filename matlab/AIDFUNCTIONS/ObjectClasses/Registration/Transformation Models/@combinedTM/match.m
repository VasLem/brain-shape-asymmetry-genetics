function out = match(obj,q,p,w)
% find the combined transform from b to a
         if isempty(obj.List),out = [];return; end
         if nargout == 1
            obj = clone(obj);
            out = obj;
         end
         if nargin<4
            w  = ones(3,TM.getNrPoints(p));
            clear tmpp;
         else
            if size(w,1)==1, w = repmat(w,3,1);end 
         end
         for i=1:1:obj.nrTM
             match(obj.List{i},q,p,w);
             inverse = match(obj.List{i},p,q,w);
             if i<obj.nrTM
                eval(inverse,q);
                q = fastClone(inverse.Evaluation);
                delete(inverse);
             end
         end
end



