    function [retS, retL, retR, retP] = ProcrustesAnova2WayAsymmetryDebuggingExperiments(X1, X2, mult)
    toCheckSizes= [100, 200, 500, 1000];
    retS = checkSize(X1, X2, toCheckSizes,mult);
    plotExp(retS, toCheckSizes, 'Population Size');
    
    %%
    toCheckLSizes= [50, 100, 150, 200, 250];
    retL = checkLSize(X1, X2, toCheckLSizes,mult);
    plotExp(retL, toCheckLSizes, 'Number of Landmarks');
    
    %%
    
    toCheckRepsNum= [1,2,3];
    retR = checkRepsNum(X1, X2, toCheckRepsNum,mult);
    plotExp(retR, toCheckRepsNum, 'Number of Replications');
    
    %%
    toCheckPermsNum= [1000, 2000, 3000];
    retP = checkPermsNum(X1, X2, toCheckPermsNum, mult);
    plotExp(retP, toCheckPermsNum, 'Number of Permutations');
    end