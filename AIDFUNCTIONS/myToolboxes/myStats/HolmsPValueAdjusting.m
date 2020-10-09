function out = HolmsPValueAdjusting(pvalues,pCrit)
% reference: Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing
         nrTest = length(pvalues);
         [sortedpvalues,sortind] = sort(pvalues,'ascend');
         good = zeros(1,nrTest);
         counter = 0;
         for i=1:1:nrTest
             counter = counter +1;
             forpCrit = pCrit/(nrTest-i+1);
             if sortedpvalues(i)<=forpCrit,good(i) = 1;end
         end
         out = sortind(good==1);
end