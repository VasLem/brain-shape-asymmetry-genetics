function [out,M] = getShapeResiduals(X,Y)
         [~,~,~,~,M] = plsregress(X,Y,min(size(X,2),size(Y,2)));
         Y_est = [ones(size(X,1),1) X]*M;
         out = Y-Y_est;
end