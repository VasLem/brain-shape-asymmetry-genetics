function varargout = refresh(obj,p,varargin)
         if nargout == 1
            obj = clone(obj);
            varargout{1} = obj;
         end
         if isempty(obj.List), return; end
         for i=1:1:obj.nrTM
             refresh(obj.List{i},p,varargin{:});
             p = obj.List{i}.Evaluation;
         end
end