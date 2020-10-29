function [obj,M] = getRIP2YModel(obj,index,Cindex)
    if nargin < 3, Cindex = [];end
    switch obj.BRIMModelType
        case 'PLSR'
            if isempty(obj.RIP2YM)
               obj.RIP2YM = nan*zeros(obj.nrX,obj.nrY);
            end
            nr2Condition = length(Cindex);
            nr2Model = length(index);
            if ~nr2Condition==0
               IndVar = [obj.RIPX(:,Cindex) obj.RIPX(:,index)];
            else
               IndVar = obj.RIPX(:,index);
            end       
            ModelInd = (nr2Condition+1:nr2Condition+nr2Model);
            DepVar = obj.Y;     
            [~,~,~,~,M] = plsregress(IndVar,DepVar,size(IndVar,2));
            obj.RIP2YM(index,:) = M(ModelInd+1,:);
        otherwise
    end
end