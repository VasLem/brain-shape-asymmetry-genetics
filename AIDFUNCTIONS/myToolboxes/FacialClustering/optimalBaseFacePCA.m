function out = optimalBaseFacePCA(betha,in,avgX,EigStd)
         
         options = psoptimset('UseParallel', false, 'CompletePoll', 'on', 'Vectorized', 'off','TolMesh',1e-4,'Display','off');   
         [out.searchmin, out.fval] = patternsearch(@(X)fiterror(X,betha,in,EigStd),avgX,[],[],[],[],[],[],options);  
         out.Y_est = [ones(size(out.searchmin,1),1) out.searchmin]*betha;   
         
end

function error = fiterror(X,M,Y,EigStd)
         Y_est = [ones(size(X,1),1) X]*M;
         %error = sum(((Y-Y_est)./EigStd).^2);
         error = sqrt(sum(((Y-Y_est)).^2));
end