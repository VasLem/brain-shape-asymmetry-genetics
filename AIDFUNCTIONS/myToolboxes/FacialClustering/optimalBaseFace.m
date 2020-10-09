function out = optimalBaseFace(betha,in,avgX)
         
         options = psoptimset('UseParallel', true, 'CompletePoll', 'on', 'Vectorized', 'off','TolMesh',5e-5,'Display','iter');   
         [out.searchmin, out.fval] = patternsearch(@(X)fiterror(X,betha,in),avgX,[],[],[],[],[],[],options);
end

function error = fiterror(X,M,Y)
         Y_est = [ones(size(X,1),1) X]*M;
         res = Y-Y_est;
         error = sqrt(mean(res.^2));
end