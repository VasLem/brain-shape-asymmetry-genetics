function out = regressionAngles(A1,B1,A2,B2,t)
tic;
[A1,B1] = eliminateNAN(A1,B1);
[A2,B2] = eliminateNAN(A2,B2);
if nargin < 5, t = 0; end
% A1, B1
    n1 = size(A1,1);
    nA = size(A1,2);
    [~,~,~,~,M1,~,~,~] = plsregress(A1,B1,nA);
% A2, B2
    n2 = size(A2,1);
    nA = size(A2,2);
    [~,~,~,~,M2,~,~,~] = plsregress(A2,B2,nA);
% Get angles 
    angles = zeros(1,nA);
    anglesd = zeros(1,nA);
    for i=1:1:nA
        [angles(i),anglesd(i)] = angle(M1(i+1,:)',M2(i+1,:)'); 
    end
    out.angles = angles;
    out.anglesd = anglesd;
if t == 0, return; end
% Random angles
%     dim = size(M1,2);
%     RandAngles = zeros(1,t);
%     parfor i=1:t
%            Dir1 = -1 + 2*rand(1,dim);Dir1 = Dir1/norm(Dir1);
%            Dir2 = -1 + 2*rand(1,dim);Dir2 = Dir2/norm(Dir2);
%            RandAngles(i) = angle(Dir1',Dir2');
%     end
%     % get P-Value 
%     RandPvalues = zeros(1,nA);
%     for i=1:1:nA
%         switch sign(angles(i));
%             case 1
%                 Count = RandAngles>=angles(i);
%             case -1
%                 Count = RandAngles<=angles(i);
%         end
%         RandPvalues(i) = sum(Count)/t;
%     end
%     out.RandPvalues = RandPvalues;
% Bootstrapping
    n = min(n1,n2);% determine the smaller group size 
    BootAngles1 = zeros(nA,t);
    BootAngles2 = zeros(nA,t);
    parfor i=1:t
       % A1 B1
       ind = randsample((1:n1),n,true);
       [~,~,~,~,M1a,~,~,~] = plsregress(A1(ind,:),B1(ind,:),nA); %#ok<*PFBNS>
       ind = randsample((1:n1),n,true);
       [~,~,~,~,M1b,~,~,~] = plsregress(A1(ind,:),B1(ind,:),nA); %#ok<*PFBNS>
       TmpAngles = zeros(nA,1);
       for k=1:1:nA
           TmpAngles(k) = angle(M1a(k+1,:)',M1b(k+1,:)'); 
       end
       BootAngles1(:,i) = TmpAngles;
       % A2, B2
       ind = randsample((1:n2),n,true);
       [~,~,~,~,M2a,~,~,~] = plsregress(A2(ind,:),B2(ind,:),nA); %#ok<*PFBNS>
       ind = randsample((1:n2),n,true);
       [~,~,~,~,M2b,~,~,~] = plsregress(A2(ind,:),B2(ind,:),nA); %#ok<*PFBNS>
       TmpAngles = zeros(nA,1);
       for k=1:1:nA
           TmpAngles(k) = angle(M2a(k+1,:)',M2b(k+1,:)'); 
       end
       BootAngles2(:,i) = TmpAngles;
    end
    for i=1:1:nA
        list = sort(BootAngles1(i,:));
        lowerci= list(0.025*t+1);
        out.
    end
    
    % get P-Value 
    BootPvalues = zeros(1,nA);
    for i=1:1:nA
        Count = abs(BootAngles(i,:))<=abs(angles(i));
        BootPvalues(i) = sum(Count)/t;
    end
    out.BootPvalues = BootPvalues;
toc;
end

function [A,B] = eliminateNAN(A,B)
         index = (1:size(A,1));
         [i,~] = find(isnan(A));
         i = unique(i);
         index = setdiff(index,i);
         A = A(index,:);
         B = B(index,:);
end

% function out = regressionAngles(A1,B1,A2,B2,t)
% tic;
% [A1,B1] = eliminateNAN(A1,B1);
% [A2,B2] = eliminateNAN(A2,B2);
% if nargin < 5, t = 0; end
% % A1, B1
%     n1 = size(A1,1);
%     nA = size(A1,2);
%     [~,~,~,~,M1,~,~,~] = plsregress(A1,B1,nA);
% % A2, B2
%     n2 = size(A2,1);
%     nA = size(A2,2);
%     [~,~,~,~,M2,~,~,~] = plsregress(A2,B2,nA);
% % Get angles 
%     angles = zeros(1,nA);
%     anglesd = zeros(1,nA);
%     for i=1:1:nA
%         [angles(i),anglesd(i)] = angle(M1(i+1,:)',M2(i+1,:)'); 
%     end
%     out.angles = angles;
%     out.anglesd = anglesd;
% if t == 0, return; end
% % Random angles
%     dim = size(M1,2);
%     RandAngles = zeros(1,t);
%     parfor i=1:t
%            Dir1 = -1 + 2*rand(1,dim);Dir1 = Dir1/norm(Dir1);
%            Dir2 = -1 + 2*rand(1,dim);Dir2 = Dir2/norm(Dir2);
%            RandAngles(i) = angle(Dir1',Dir2');
%     end
%     % get P-Value 
%     RandPvalues = zeros(1,nA);
%     for i=1:1:nA
%         switch sign(angles(i));
%             case 1
%                 Count = RandAngles>=angles(i);
%             case -1
%                 Count = RandAngles<=angles(i);
%         end
%         RandPvalues(i) = sum(Count)/t;
%     end
%     out.RandPvalues = RandPvalues;
% % Bootstrapping
%     BootAngles = zeros(nA,t);
%     parfor i=1:t
%        % A1, B1
%        ind = randsample((1:n1),n1,true);
%        [~,~,~,~,M1for,~,~,~] = plsregress(A1(ind,:),B1(ind,:),nA); %#ok<*PFBNS>
%        % A2, B2
%        ind = randsample((1:n2),n2,true);
%        [~,~,~,~,M2for,~,~,~] = plsregress(A2(ind,:),B2(ind,:),nA);
%        TmpAngles = zeros(nA,1);
%        for k=1:1:nA
%            TmpAngles(k) = angle(M1for(k+1,:)',M2for(k+1,:)'); 
%        end
%        BootAngles(:,i) = TmpAngles;
%     end
%     % get P-Value 
%     BootPvalues = zeros(1,nA);
%     for i=1:1:nA
%         Count = abs(BootAngles(i,:))<=abs(angles(i));
%         BootPvalues(i) = sum(Count)/t;
%     end
%     out.BootPvalues = BootPvalues;
% toc;
% end
% 
% function [A,B] = eliminateNAN(A,B)
%          index = (1:size(A,1));
%          [i,~] = find(isnan(A));
%          i = unique(i);
%          index = setdiff(index,i);
%          A = A(index,:);
%          B = B(index,:);
% end