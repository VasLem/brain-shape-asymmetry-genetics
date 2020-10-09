function xoverKids  = crossoverheuristic(obj)
%CROSSOVERHEURISTIC Move from worst parent to slightly past best parent.
%   XOVERKIDS = CROSSOVERHEURISTIC(PARENTS,OPTIONS,GENOMELENGTH, ...
%   FITNESSFCN,THISSCORE,THISPOPULATION,RATIO) creates the crossover 
%   children XOVERKIDS of the given population THISPOPULATION using the 
%   available PARENTS, the current scores THISSCORE and a RATIO. This kind
%   of recombination will only work for real valued genomes.
%
%   Example:
%    Create an options structure using CROSSOVERHEURISTIC as the crossover
%    function with default ratio of 1.2
%     options = gaoptimset('CrossoverFcn' , @crossoverheuristic);
%
%    Create an options structure using CROSSOVERHEURISTIC as the crossover
%    function with RATIO of 1.5
%     ratio = 1.5;
%     options = gaoptimset('CrossoverFcn' , {@crossoverheuristic,ratio});

%   Copyright 2003-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/08/29 08:24:18 $

parents = obj.CrossoverParents;

% here's the default if nothing is passed in
if nargin < 7 || isempty(obj.CrossoverHeuristicRatio)
    obj.CrossoverHeuristicRatio = 1.2;
end
% How many children am I being asked to produce?
nKids = length(parents)/2;
% Allocate space for the kids
xoverKids = zeros(nKids,obj.Dim);
% To move through the parents twice as fast as thekids are
% being produced, a separate index for the parents is needed
index = 1;
for i=1:nKids
    % get parents
    parent1 = obj.Population(parents(index),:);
    score1 = obj.Scores(parents(index));
    index = index + 1;
    parent2 = obj.Population(parents(index),:);
    score2 = obj.Scores(parents(index));
    index = index + 1;
    
    % move a little past the better away from the worst
    if(score1 < score2) % parent1 is the better of the pair
        xoverKids(i,:) = parent2 + obj.CrossoverHeuristicRatio .* (parent1 - parent2);
    else % parent2 is the better one
        xoverKids(i,:) = parent1 + obj.CrossoverHeuristicRatio .* (parent2 - parent1);
    end
    % Make sure that offspring are feasible w.r.t. linear constraints
    %xoverKids(i,:) = trimKid(obj,xoverKids(i,:));
    if ~feasibleKid(obj,xoverKids(i,:))
       % Children are arithmetic mean of two parents (feasible w.r.t
       % linear constraints)
       alpha = rand;
       xoverKids(i,:) = alpha*obj.Population(r1,:) + ...
                (1-alpha)*obj.Population(r2,:); 
    end
end
