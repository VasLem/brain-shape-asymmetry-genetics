function out = derivate(obj,p)
      if ~isField(obj,obj.ActiveField)||isempty(obj.List)
         if nargout == 1, out = [];end
         return;
      end
      % current evaluation is different for sure than what is needed
      % if nrV are not equal and/or obj.Evaluation is invalid, seen in
      % obj.nrV
      % Adviced to: always eval before derivate in main procedures!!!
      % Calling eval always, here is computational un-efficient!!!!!!
      nrV = TM.getNrPoints(p);
      if ~(nrV==obj.nrV), eval(obj,p); end 
      % current ActiveModel, function of ActiveP
      % sets the ActiveP of current ActiveModel
      mindex = obj.ActiveModel;
      % if not the last model in the list, take the evaluation of the
      % next model in the list to start with
      if mindex < obj.nrTM
         p = obj.List{mindex+1}.Evaluation;
      end
      % perform derivate
      grad = derivate(obj.List{mindex},p);
      % can be empty when derivating towards a parameter not having an
      % effect on the current ActiveField
      if ~isempty(grad)
          % if not the first model in the list, ripple the grad back to front
          % using the appropriate chainrules
          if mindex>1
             for k=mindex-1:-1:1
                 grad = chainRule(obj.List{k},grad,obj.List{k+1}.Evaluation);
             end
          end
      end
      if nargout == 1, out = grad; return; end
      % Storing approriatly based on ActiveField see set.Derivative in
      % TM class
      obj.Derivative = grad;
end


% function [mindex,pindex] = geListIndex(index,nrP,nrM)
%          nr = 0;
%          for mindex = 1:1:nrM
%              nr = nr + nrP(mindex);
%              if index<nr
%                  break;
%              end
%          end
%          pindex = nrP(mindex)-(nr-index);         
% end

%       if isempty(obj.List), return; end
%       nr = size(p,2);
%       nrP = obj.nrTMP;
%       if ~(nr==obj.nrE),clear(obj);eval(obj,p,varargin{:});end
%       Input = find(strcmp(varargin,'PIndex'));
%       if ~isempty(Input)% memory efficient
%          [mindex,pindex] = geListIndex(varargin{Input+1},nrP,obj.nrTM);
%          if mindex < obj.nrTM
%                 p = obj.List{mindex+1}.Evaluation;
%          end
%          grad = derivate(obj.List{mindex},p,'Pindex',pindex,varargin{:});
%          if mindex>1
%             for k=mindex-1:-1:1
%                 grad = chainRule(obj.List{k},grad,varargin{:});
%             end
%          end
%       else % computation efficient
%           grad = zeros(3,nr,obj.nrP);
%           index = obj.nrP;
%           for i=obj.nrTM:-1:1
%               tmpgrad = derivate(obj.List,p,varargin{:});
%               if ~(i==1)
%                  tmpgrad = matrix2list(tmpgrad,3); 
%                  for k=i-1:-1:1
%                      tmpgrad = chainRule(obj.List{k},tmpgrad,varargin{:});
%                  end
%                  tmpgrad = list2matrix(tmpgrad,3,nrP(i));
%               end
%               grad(:,:,index:index-nrP(i)+1) = tmpgrad;
%               p = obj.List{i}.Evaluation;
%               index = index-nrP(i);         
%           end
%       end
%       if nargout == 1, varargout{1} = grad; return; end
%       obj.Grads = grad;
