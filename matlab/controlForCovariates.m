%%TODO
function [controlled,  expVariance]  = controlForCovariates(covariates, toBeControlled)
    [~,~,~,~,M, PCTVAR] = plsregress(covariates,toBeControlled,min(size(covariates,2),size(toBeControlled,2)));
    expVariance = sum(PCTVAR, 2);
    expVariance = expVariance(2);
    Y_est = [ones(size(covariates,1),1) covariates]*M;
    controlled =toBeControlled - Y_est;
end