function  [I,D,F] = AMMIAnalysis(X1, X2)
    X = [X1(:) X2(:)];
    I = X - mean(X,1);
    D = X - mean(X,2);
    [F,score,latent,tsquared,explained,mu] = pca(X - mean(X,1) - mean(X,2));
    
    
    

end