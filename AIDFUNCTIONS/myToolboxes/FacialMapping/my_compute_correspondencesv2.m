function [correspondingFeatures,correspondingFlags,pullF] = my_compute_correspondencesv2(floatingFeatures, targetFeatures,floatingFlags,targetFlags,correspondencesSymmetric,correspondencesNumNeighbours,correspondencesFlagThreshold)
         % extracting points from features
             floatingPoints = floatingFeatures(:,1:3);
             nF = size(floatingPoints,1);
             targetPoints = targetFeatures(:,1:3);
             nT = size(targetPoints,1);
             %correspondingFeatures = floatingFeatures;
             nFeatures = size(floatingFeatures,2);
         % PUSH corresponding Points
             [pushK,pushD] = knnsearch(targetFeatures,floatingFeatures,'K',correspondencesNumNeighbours);
             pushW = 1./(pushD.^2);  
             pushA = zeros(nF,nT);
             for k=1:correspondencesNumNeighbours
                 ind = sub2ind([nF nT],1:nF,pushK(:,k)');
                 pushA(ind) = pushW(:,k);
             end
             sumPushA = sum(pushA,2);
             %pushC = pushA*targetPoints./repmat(sumPushA,1,3);   
             pushF = pushA*targetFlags./sumPushA;           
             if ~correspondencesSymmetric            
                 %correspondingFeatures(:,1:3) = pushA*targetPoints./repmat(sumPushA,1,3);
                 correspondingFeatures = pushA*targetFeatures./repmat(sumPushA,1,nFeatures);
                 correspondingFlags = single(pushF>=correspondencesFlagThreshold);
                 return;
             end   
         % PULL corresponding Points
             [pullK,pullD] = knnsearch(floatingFeatures,targetFeatures,'K',correspondencesNumNeighbours);
             pullW = 1./(pullD.^2);
             pullA = zeros(nT,nF);
             for k=1:correspondencesNumNeighbours
                 ind = sub2ind([nT nF],1:nT,pullK(:,k)');
                 pullA(ind) = pullW(:,k);
             end
             sumPullA = sum(pullA,2);
%              pullC = pullA*floatingPoints./repmat(sumPullA,1,3);
             pullF = pullA*floatingFlags./sumPullA;
         % RELATING PUSH AND PULL MESH % Makes things slower but is more
%          % correct 
%              [rK,rD] = knnsearch(pullC,floatingFeatures(:,1:3),'K',correspondencesNumNeighbours);
%              rW = 1./(rD.^2);
%              rA = zeros(nF,nT);
%              for k=1:correspondencesNumNeighbours
%                  ind = sub2ind([nF nT],1:nF,rK(:,k)');
%                  rA(ind) = rW(:,k);
%              end
%              sumrA = sum(rA,2);
%              pullF3 = rA*pullF./sumrA;
%              %pullF2 = pushA*pullF./sumPushA;
         % COMBINING both
         
         pushA = pushA./repmat(sum(pushA,2),1,nT);
         pullA = pullA./repmat(sumPullA,1,nF); 
%          
         symA = pushA+pullA';
         %symC = symA*targetPoints./repmat(sum(symA,2),1,3);
         
%          symF = pullF3.*pushF;
         
         pullF2 = symA*pullF./sum(symA,2);
         
         symF = pullF2.*pushF;
         
         %correspondingFeatures = floatingFeatures;
         %correspondingFeatures(:,1:3) = symC;
         correspondingFeatures = symA*targetFeatures./repmat(sum(symA,2),1,nFeatures);
         correspondingFlags = single(symF>=correspondencesFlagThreshold);
         
         
         
         
                 
end


% v = viewer(floatingMesh);
% floatingMesh.Value = symF;
% floatingMesh.ColorMode = 'Indexed';
% symMesh.Axes = v.RenderAxes;
% symMesh.Visible = true;
% symMesh.SingleColor = [0.8 0.8 0.8];
% 
% diff = vDifferences(symMesh,floatingMesh);
% dist = vDistances(symMesh,floatingMesh);
% quiver3(v.RenderAxes,floatingMesh.Vertices(1,:)',floatingMesh.Vertices(2,:)',floatingMesh.Vertices(3,:)',...
% diff(1,:)',diff(2,:)',diff(3,:)',6,'color',[1 1 1]);


% function [correspondingFeatures,correspondingFlags] = my_compute_correspondences(floatingFeatures, targetFeatures,floatingFlags,targetFlags,correspondencesSymmetric,correspondencesNumNeighbours,correspondencesFlagThreshold)
%          %% now with affinity matrices
%          
%          % settings 
%          correspondencesSymmetric = true;
%          correspondencesNumNeighbours = 5;
%          correspondencesFlagThreshold = 0.9;
%          
%          pushMesh = clone(floatingMesh);
%          pullMesh = clone(targetMesh);
%          symMesh = clone(floatingMesh);
%          
%          % extracting points from features
%          floatingPoints = floatingFeatures(:,1:3);
%          nF = size(floatingPoints,1);
%          targetPoints = targetFeatures(:,1:3);
%          nT = size(targetPoints,1);
%          % PUSH corresponding Points
%          [pushK,pushD] = knnsearch(targetFeatures,floatingFeatures,'K',correspondencesNumNeighbours);
%          pushW = 1./(pushD.^2);
%          
%          pushA = zeros(nF,nT);
%          for k=1:correspondencesNumNeighbours
%              ind = sub2ind([nF nT],1:nF,pushK(:,k)');
%              pushA(ind) = pushW(:,k);
%          end
%          sumPushA = sum(pushA,2);
%          pushC = pushA*targetPoints./repmat(sumPushA,1,3);   
%          pushF = pushA*targetFlags./sumPushA;
%          
%          pushMesh.Vertices = pushC';
%          pushMesh.Value = pushF;
%          viewer(pushMesh);
%          
%          % PULL corresponding Points
%          [pullK,pullD] = knnsearch(floatingFeatures,targetFeatures,'K',correspondencesNumNeighbours);
%          pullW = 1./(pullD.^2);
%          
%          pullA = zeros(nT,nF);
%          for k=1:correspondencesNumNeighbours
%              ind = sub2ind([nT nF],1:nT,pullK(:,k)');
%              pullA(ind) = pullW(:,k);
%          end
%          sumPullA = sum(pullA,2);
%          pullC = pullA*floatingPoints./repmat(sumPullA,1,3);
%          pullF = pullA*floatingFlags./sumPullA;
%          pullMesh.Vertices = pullC';
%          pullMesh.Value = pullF;
%          viewer(pullMesh);
%          
%          % RELATING PUSH AND PULL MESH % Makes things slower but is more
%          % correct 
%          [rK,rD] = knnsearch(pullC,floatingFeatures(:,1:3),'K',correspondencesNumNeighbours);
%          rW = 1./(rD.^2);
%          rA = zeros(nF,nT);
%          for k=1:correspondencesNumNeighbours
%              ind = sub2ind([nF nT],1:nF,rK(:,k)');
%              rA(ind) = rW(:,k);
%          end
%          sumrA = sum(rA,2);
%          pullF3 = rA*pullF./sumrA;
%          pullF2 = pushA*pullF./sumPushA;
%          %% COMBINING both
%          
%          %pushA = pushA./repmat(sum(pushA,2),1,nT);
%          %pullA = pullA./repmat(sumPullA,1,nF); 
%          
%          symA = pushA+pullA';
%          symC = symA*targetPoints./repmat(sum(symA,2),1,3);        
%          
%          symMesh.Vertices = symC';
%          viewer(symMesh);
%          
%          symF = pullF3.*pushF;
%          %symF = symA*pullF./sum(symA,2).*pushF;
%          %symMesh.Value = double(symF>0.99)';
%          symMesh.Value = symF;        
% end
% 
% 
% v = viewer(floatingMesh);
% floatingMesh.Value = symF;
% floatingMesh.ColorMode = 'Indexed';
% symMesh.Axes = v.RenderAxes;
% symMesh.Visible = true;
% symMesh.SingleColor = [0.8 0.8 0.8];
% 
% diff = vDifferences(symMesh,floatingMesh);
% dist = vDistances(symMesh,floatingMesh);
% quiver3(v.RenderAxes,floatingMesh.Vertices(1,:)',floatingMesh.Vertices(2,:)',floatingMesh.Vertices(3,:)',...
% diff(1,:)',diff(2,:)',diff(3,:)',6,'color',[1 1 1]);



% 
%          % settings 
%          correspondencesSymmetric = true;
%          correspondencesNumNeighbours = 5;
%          correspondencesFlagThreshold = 0.9;
%          
%          pushMesh = clone(floatingMesh);
%          pullMesh = clone(targetMesh);
%          pullForwardMesh = clone(floatingMesh);
%          
%          % extracting points from features
%          floatingPoints = floatingFeatures(:,1:3);
%          nF = size(floatingPoints,1);
%          targetPoints = targetFeatures(:,1:3);
%          nT = size(targetPoints,1);
%          % PUSH corresponding Points
%          [pushK,pushD] = knnsearch(targetFeatures,floatingFeatures,'K',correspondencesNumNeighbours);
%          pushW = 1./(pushD.^2);
%          
%          pushC = floatingPoints;
%          
%          parfor i=1:1:nF
%              %i=1;
%              tmpP = targetPoints(pushK(i,:),:);
%              tmpW = repmat(pushW(i,:)',1,3);
%              pushC(i,:) = sum(tmpP.*tmpW,1)/sum(pushW(i,:));  
%          end
%          pushMesh.Vertices = pushC';
%          
%          if ~correspondencesSymmetric
%             correspondingFeatures(:,1:3) = pushC;
%             % to do flag handling
%             return;
%          end
%              
%          % PULL corresponding Points
%          [pullK,pullD] = knnsearch(floatingFeatures,targetFeatures,'K',correspondencesNumNeighbours);
%          pullW = 1./(pullD.^2);
%          pullC = targetPoints;      
%          parfor i=1:1:nT
%              %disp(num2str(i));
%              tmpP = floatingPoints(pullK(i,:),:);
%              tmpW = repmat(pullW(i,:)',1,3);
%              pullC(i,:) = sum(tmpP.*tmpW,1)/sum(pullW(i,:));  
%          end
%          pullMesh.Vertices = pullC';
%          
%          
%          pullForwardC = floatingPoints;
%          
%          for d=1:1:3
%              %d = 1;
%              F = scatteredInterpolant(double(pullC),double(targetPoints(:,d)),'nearest');
%              pullForwardC(:,d) = F(double(floatingPoints)); 
%          end
%          pullForwardMesh.Vertices = pullForwardC';
%          
%          finalC = (pushC+pullForwardC)/2;
%          finalMesh = clone(floatingMesh);
%          finalMesh.Vertices = finalC';
         
         


