function out = derivate(obj,p)
      if ~isField(obj,obj.ActiveField)
         if nargout == 1, out = [];end
         return;
      end
      p = TM.getPoints(p);
      nr = size(p,2);
      switch obj.ActiveP
          case 1              
              if isempty(obj.c),obj.c = p;end
              p = p - repmat(obj.c,1,nr);
              grad = dTdRx(obj)*p;
          case 2
              if isempty(obj.c),obj.c = p;end
              p = p - repmat(obj.c,1,nr);
              grad = dTdRy(obj)*p;
          case 3
              if isempty(obj.c),obj.c = p;end
              p = p - repmat(obj.c,1,nr);
              grad = dTdRz(obj)*p;
          case 4
              grad = repmat(dTdTx(obj),1,nr);
          case 5
              grad = repmat(dTdTy(obj),1,nr);
          case 6
              grad = repmat(dTdTz(obj),1,nr);
          otherwise
              if nargout == 1, out = [];end
              return;
      end
      if nargout==1, out = grad;return; end
      obj.Derivative = grad;
end 
%       p = getP(p);
%       if ~size(p,1)==3, p = p'; end
%       nr = size(p,2);
%       if isempty(obj.c),obj.c = p;end
%       p = p - repmat(obj.c,1,nr);
%       Input = find(strcmp(varargin,'PIndex'));
%       if ~isempty(Input)% memory efficient 
%           switch varargin{Input+1}
%               case 1
%                   grad = dTdRx(obj)*p;
%               case 2
%                   grad = dTdRy(obj)*p;
%               case 3
%                   grad = dTdRz(obj)*p;
%               case 4
%                   grad = repmat(dTdTx(obj),1,nr);
%               case 5
%                   grad = repmat(dTdTy(obj),1,nr);
%               case 6
%                   grad = repmat(dTdTz(obj),1,nr);
%               otherwise
%           end
%       else  % computation efficient    
%           grad = zeros(3,nr,obj.nrP);
%           grad(:,:,1) = dTdRx(obj)*p;
%           grad(:,:,2) = dTdRy(obj)*p;
%           grad(:,:,3) = dTdRz(obj)*p;        
%           grad(:,:,4) = repmat(dTdTx(obj),1,nr);
%           grad(:,:,5) = repmat(dTdTy(obj),1,nr);
%           grad(:,:,6) = repmat(dTdTz(obj),1,nr);
%       end
%       if nargout==1, varargout{1} = grad;return; end
%       obj.Derivative = grad;
