function [Size,Shape,Space] = myGetBrainClusterSizeShape(data,adjustforsize,percvar,covariates)
    if nargin< 2, adjustforsize = true; end
    if nargin<3, percvar = 90; end
    if nargin <4, covariates = []; end
    [nObs,nVoxels] = size(data);
    if ~isempty(covariates)
       avgD = mean(data,1); 
       res = getResiduals(covariates,data);
       data = repmat(avgD,nObs,1)+res;
    end
    Size = sum(data,2)./nVoxels;
    if adjustforsize == true
        data = data./repmat(Size,1,nVoxels);
    end
    Space = plainPCA;
    getAverage(Space,data');
    getModel(Space,data');
    stripPercVar(Space,percvar);
    Shape = Space.Tcoeff./repmat(Space.EigStd',Space.n,1);   
end

% 
% Space = plainPCA;
% getAverage(Space,res');
% getModel(Space,res');
% new = stripPercVar(Space,percvar);