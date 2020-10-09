function [LocalE,LocalS,GlobalS] = getRegressionStatistic(M,Res,n,nA,nB,type)
         switch lower(type)
             case 'value'
                 LocalE = zeros(nA,nB);
                 LocalS = zeros(nA,nB);
                 GlobalS = zeros(nA,1);
                 for i=1:1:nA
                     LocalE(i,:) = M(1+i,:);
                     LocalS(i,:) = abs(LocalE(i,:));
                     GlobalS(i) = sum(LocalS(i,:));
                 end
             case 'shape'
                 LocalE = zeros(3,nB/3,nA);
                 LocalS = zeros(nA,nB/3);
                 GlobalS = zeros(nA,1);
                 for i=1:1:nA
                     LocalE(:,:,i) = reshape(M(1+i,:),3,nB/3);
                     SS = sum(LocalE(:,:,i).^2);
                     LocalS(i,:) = sqrt(SS);
                     GlobalS(i) = sum(SS);
                 end 
             otherwise
                 error('unknown type');
         end        

end

% function [LocalE,LocalS,GlobalS] = getRegressionStatistic(M,Res,nA,nB,type)
%          switch lower(type)
%              case 'value'
%                  LocalE = zeros(nA,nB);
%                  LocalS = zeros(nA,nB);
%                  GlobalS = zeros(nA,1);
%                  for i=1:1:nA
%                      LocalE(i,:) = M(1+i,:);
%                      LocalS(i,:) = abs(LocalE(i,:));
%                      GlobalS(i) = sum(LocalS(i,:));
%                  end
%              case 'shape'
%                  LocalE = zeros(3,nB/3,nA);
%                  LocalS = zeros(nA,nB/3);
%                  GlobalS = zeros(nA,1);
%                  for i=1:1:nA
%                      LocalE(:,:,i) = reshape(M(1+i,:),3,nB/3);
%                      SS = sum(LocalE(:,:,i).^2);
%                      LocalS(i,:) = sqrt(SS);
%                      GlobalS(i) = sum(SS);
%                  end 
%              otherwise
%                  error('unknown type');
%          end        
% 
% end
