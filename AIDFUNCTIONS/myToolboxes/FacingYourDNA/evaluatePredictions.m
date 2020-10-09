function out = evaluatePredictions(TestY,Y_est)
         [nS,dim] = size(TestY);
%          ANGLES = nan*zeros(nS,nS);
%          for i=1:1:nS
%              %i=1;
%              for j=1:1:nS
%                  %j=1;
%                  ANGLES(i,j) = angle(TestY(i,:)',Y_est(j,:)');
%              end
%          end
         
         
%          TrueAngles = diag(ANGLES);
%          FalseAngles = triu(ANGLES);
%          
         DIST = pdist2(Y_est,TestY,'euclidean');
           
         R = zeros(1,nS);
         for i=1:1:nS
            distances = DIST(i,:);
            tdistance = distances(i);
            fdistance = distances(setdiff(1:nS,i));
            R(i) = sum(fdistance<=tdistance);
         end
         
         R = (R./nS).*100;
         CumRanks = zeros(1,100);
         for cr=1:1:100
             CumRanks(:,cr) = (sum(R<=cr)/nS).*100;
         end
         
         figure;set(gca,'ylim',[0 100],'xlim',[1 100]);hold on;
         xlabel('Rank (%)');ylabel('Identified (%)' );grid on;
         title('Identification');
         plot(1:100,1:100,'k--');
         plot(1:100,squeeze(CumRanks(1,:)),'b-','LineWidth',2);
         
         DIST = pdist2(TestY,Y_est,'cosine');
         
         R = zeros(1,nS);
         for i=1:1:nS
            distances = DIST(i,:);
            tdistance = distances(i);
            fdistance = distances(setdiff(1:nS,i));
            R(i) = sum(fdistance<=tdistance);
         end
         
         R = (R./nS).*100;
         CumRanks = zeros(1,100);
         for cr=1:1:100
             CumRanks(:,cr) = (sum(R<=cr)/nS).*100;
         end
         plot(1:100,squeeze(CumRanks(1,:)),'g-','LineWidth',2);
         
         out = [];

end