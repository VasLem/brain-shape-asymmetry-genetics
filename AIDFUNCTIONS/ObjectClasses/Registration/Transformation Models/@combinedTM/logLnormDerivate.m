function out = logLnormDerivate(obj)
     if isempty(obj.List), return; end
     %obj.ActiveModel, set the correct ActiveP for the relevant List
     lnormgrad = logLnormDerivate(obj.List{obj.ActiveModel});
     if nargout == 1, out = lnormgrad; return; end
     obj.LnormDerivative = lnormgrad;    
end
     

% nr = obj.nrTMP;
%      Input = find(strcmp(varargin,'PIndex'));
%      if ~isempty(Input)% memory efficient
%         [mindex,pindex] = geListIndex(varargin{Input+1},nr,obj.nrTM);
%         lnormgrad = logLnormDerivate(obj.List(mindex),'Pindex',pindex,varargin{:});
%      else
%         lnormgrad = zeros(size(obj.nrP));
%         index = 1;
%         for i=1:1:obj.nrTM
%             lnormgrad(index:index+nr(i)-1) = logLnormDerivate(obj.List{i},varargin{:});
%             index = nr(i)+1;       
%         end
%      end
%      if nargout == 1, varargout{1} = lnormgrad;return; end
%      obj.LnormDerivative = lnormgrad;