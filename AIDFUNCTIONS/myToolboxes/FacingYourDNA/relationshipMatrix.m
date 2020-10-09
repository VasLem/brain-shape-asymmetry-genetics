function GD = relationshipMatrix(in)
         [nrS,nrSNP] = size(in);
         
%          freq = zeros(1,nrSNP);
%          for i=1:1:nrSNP
%              freq(i) = 
%              
%          end
         
         
         GD = zeros(nrS,nrS);
         for i=1:1:nrS
            for j=i+1:1:nrS
                GD(i,j)= getASDistance(in(i,:),in(j,:));
            end
         end
         



end

