% This function can be used to aggregate pvalues
% retrieved from:
% https://imaging.mrc-cbu.cam.ac.uk/statswiki/FAQ/CombiningPvalues

function pcomb = stouffer(p)
% Stouffer et al's (1949) unweighted method for combination of 
% independent p-values via z's 
	pcomb = (erfc(sum(sqrt(2) * erfcinv(2.*p),2)/sqrt(2*size(p,2))))/2;    
end
