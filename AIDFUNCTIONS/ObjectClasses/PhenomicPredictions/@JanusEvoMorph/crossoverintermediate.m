function xoverKids  = crossoverintermediate(obj)
%CROSSOVERINTERMEDIATE Weighted average of the parents.
%   XOVERKIDS = CROSSOVERINTERMEDIATE(PARENTS,OPTIONS,GENOMELENGTH, ...
%   FITNESSFCN,UNUSED,THISPOPULATION,RATIO) creates the crossover children
%   XOVERKIDS of the given population THISPOPULATION using the available
%   PARENTS. This kind of recombination can only work for real valued
%   genomes. A point is chosen partway between the two parents.
%
%   RATIO controls the range over which the children can occur.
%   A scalar ratio will generate children on the line between the parents.
%   A vector ratio (length = GenomeLength) can produce children anywhere in
%   the hypercube with the parents at opposite corners.
%
%   Example:
%    Create an options structure using CROSSOVERINTERMEDIATE as the crossover
%    function with default ratio = ones(1,GenomeLength)
%     options = gaoptimset('CrossoverFcn' ,@crossoverintermediate);
%   Create an options structure using CROSSOVERINTERMEDIATE as the crossover
%    function with RATIO of .5
%     ratio = 0.5
%     options = gaoptimset('CrossoverFcn' ,{@crossoverintermediate, ratio});
%

%   Copyright 2003-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/08/29 08:24:19 $
parents = obj.CrossoverParents;
% here's the default if nothing is passed in.
if nargin < 7 || isempty(obj.CrossoverRatio)
    ratio = ones(1,obj.Dim);
end

% How many children to produce?
nKids = length(parents)/2;
% Allocate space for the kids
xoverKids = zeros(nKids,obj.Dim);

% To move through the parents twice as fast as thekids are
% being produced, a separate index for the parents is needed
index = 1;

for i=1:nKids
    % get parents
    parent1 = obj.Population(parents(index),:);
    index = index + 1;
    parent2 = obj.Population(parents(index),:);
    index = index + 1;

    % a random number (or vector) on the interval [0,ratio]
    scale = obj.CrossoverIntermediateRatio .* rand(1,length(obj.CrossoverIntermediateRatio));

    % a child gets half from each parent
    xoverKids(i,:) = parent1 + scale .* (parent2 - parent1);
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