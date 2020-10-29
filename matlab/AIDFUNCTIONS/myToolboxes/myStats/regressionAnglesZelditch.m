function out = regressionAnglesZelditch(A1,B1,A2,B2,t)
% this implementation is based on the book of Zelditch 2004, p252
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
    dim = size(M1,2);
    RandAngles = zeros(t,nA,2);
    parfor i=1:t
           % First distribution of random angles
           tmpAngles = zeros(nA,2);
           for k=1:1:nA
            tmp = M1(k,:);
            ind = randperm(dim);Dir1 = tmp(ind);
            ind = randperm(dim);Dir2 = tmp(ind);
            tmpAngles(k,1) = angle(Dir1',Dir2');
           end
           for k=1:1:nA
            tmp = M2(k,:);
            ind = randperm(dim);Dir1 = tmp(ind);
            ind = randperm(dim);Dir2 = tmp(ind);
            tmpAngles(k,2) = angle(Dir1',Dir2');
           end
           RandAngles(i,:,:) = tmpAngles;
    end
    % get P-Value 
    RandPvalues = zeros(1,nA);
    for i=1:1:nA
        switch sign(angles(i));
            case 1
                Count = [RandAngles(:,i,1); RandAngles(:,i,2)] >=angles(i);
            case -1
                Count = [RandAngles(:,i,1); RandAngles(:,i,2)] <=angles(i);
        end
        RandPvalues(i) = sum(Count)/(2*t);
    end
    out.RandPvalues = RandPvalues;
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
    out.AnglesBoot = zeros(size(out.angles));
    for i=1:1:nA
        tmp = 0;
        list = sort((acosd(BootAngles1(i,:))));
        upperci1(i) = list(0.975*t);
        if (out.anglesd(i))>upperci1(i), tmp = tmp + 1; end
        list = sort((acosd(BootAngles2(i,:))));
        upperci2(i) = list(0.975*t);
        if (out.anglesd(i))>upperci2(i), tmp = tmp + 1; end
        out.AnglesBoot(i) = tmp;
    end
    out.Boot1 = upperci1;
    out.Boot2 = upperci2;
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

function [out,outd] = angle(v1,v2)
      T = v1'*v2;
      N = sqrt((v1'*v1)*(v2'*v2));
      out = T/N;
      outd = acosd(out);
end