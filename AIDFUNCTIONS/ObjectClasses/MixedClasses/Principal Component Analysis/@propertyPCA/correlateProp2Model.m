function [R,P,S,pc] = correlateProp2Model(obj,propvalues,s)
         if isempty(obj.Model), R = [];return; end
         if obj.nrP==0; R = []; return; end
         coeff = obj.Model.Tcoeff;
         if nargin < 3, s = 0.05; end
         data = [coeff propvalues'];
         [R,P] = corrcoef(data);
         R = R(1:obj.nrMC,end-length(prop)+1:end);
         P = P(1:obj.nrMC,end-length(prop)+1:end);
         S = zeros(size(R));
         S(find(P<s)) = 1; %#ok<FNDSB>
         pc = find(S);
end