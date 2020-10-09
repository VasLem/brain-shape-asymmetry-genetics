function out = FalsDiscoveryRateControlling(pvalues,pCrit)
% reference: Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing

%          nrTest = length(pvalues);
%          [sortedpvalues,sortind] = sort(pvalues,'ascend');
%          counter = 0;
%          for i=1:1:nrTest
%              counter = counter +1;
%              forpCrit = (i/nrTest)*pCrit;
%              if sortedpvalues(i)>forpCrit; break;end
%          end
%          out = sortind(1:i-1);
%          
         
         nrTest = length(pvalues);
         [sortedpvalues,sortind] = sort(pvalues,'descend');
         counter = 0;
         for i=1:1:nrTest
             counter = counter +1;
             forpCrit = (i/nrTest)*pCrit;
             if sortedpvalues(i)<=forpCrit;break;end
         end
         out = sortind(i:end);
         
         
         
         
%          nrTest = length(pvalues);
%          [sortedpvalues,sortind] = sort(pvalues,'descend');
%          counter = 0;
%          for i=1:1:nrTest
%              counter = counter +1;
%              forpCrit = ((nrTest-i+1)/nrTest)*pCrit;
%              if sortedpvalues(i)<=forpCrit; break;end
%          end
%          out = sortind(counter:end);
end