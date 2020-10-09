function xoverKids  = crossoverscattered(obj)
%CROSSOVERSCATTERED Position independent crossover function.
%   XOVERKIDS = CROSSOVERSCATTERED(PARENTS,OPTIONS,GENOMELENGTH, ...
%   FITNESSFCN,SCORES,THISPOPULATION) creates the crossover children XOVERKIDS
%   of the given population THISPOPULATION using the available PARENTS.
%   In single or double point crossover, genomes that are near each other tend
%   to survive together, whereas genomes that are far apart tend to be
%   separated. The technique used here eliminates that effect. Each gene has an
%   equal chance of coming from either parent. This is sometimes called uniform
%   or random crossover.
%
%   Example:
%    Create an options structure using CROSSOVERSCATTERED as the crossover
%    function
%     options = gaoptimset('CrossoverFcn' ,@crossoverscattered);

%   Copyright 2003-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/08/29 08:24:20 $
    parents = obj.CrossoverParents;
    % How many children to produce?
    nKids = length(parents)/2;
    % Allocate space for the kids
    xoverKids = zeros(nKids,obj.Dim);
    % for each kid...
    P1 = parents(1:nKids);
    P1 = obj.Population(P1,:);
    P2 = parents(nKids+1:end);
    P2 = obj.Population(P2,:);
    test = rand(nKids,obj.Dim)>0.5;
    xoverKids(test) = P1(test);
    xoverKids(~test) = P2(~test);
end

%% MEMORY INEFFICIENT IMPLEMENTATION
% for i=1:nKids
%     % get parents
%     r1 = parents(index);
%     index = index + 1;
%     r2 = parents(index);
%     index = index + 1;
%     % Randomly select half of the genes from each parent
%     % This loop may seem like brute force, but it is twice as fast as the
%     % vectorized version, because it does no allocation.
%     for j = 1:obj.Dim
%         if(rand > 0.5)
%             xoverKids(i,j) = obj.Population(r1,j);
%         else
%             xoverKids(i,j) = obj.Population(r2,j);
%         end
%     end
%     % Make sure that offspring are feasible w.r.t. linear constraints
%     %xoverKids(i,:) = trimKid(obj,xoverKids(i,:));
%     if ~feasibleKid(obj,xoverKids(i,:))
%        % Children are arithmetic mean of two parents (feasible w.r.t
%        % linear constraints)
%        alpha = rand;
%        xoverKids(i,:) = alpha*obj.Population(r1,:) + ...
%                 (1-alpha)*obj.Population(r2,:); 
%     end
% end
