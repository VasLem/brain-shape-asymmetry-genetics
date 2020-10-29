function mutationChildren = mutationuniform(obj)
%MUTATIONUNIFORM Uniform multi-point mutation.
%   MUTATIONCHILDREN = MUTATIONUNIFORM(PARENTS,OPTIONS,GENOMELENGTH,...
%                      FITNESSFCN,STATE,THISSCORE,THISPOPULATION, ...
%                      MUTATIONRATE) Creates the mutated children using
%   uniform mutations at multiple points. Mutated genes are uniformly 
%   distributed over the range of the gene. The new value is NOT a function
%   of the parents value for the gene.
%
%   Example:
%     options = gaoptimset('MutationFcn', @mutationuniform); 
%
%   This will create an options structure specifying the mutation
%   function to be used is MUTATIONUNIFORM.  Since the MUTATIONRATE is
%   not specified, the default value of 0.01 is used.
%
%     mutationRate = 0.05;
%     options = gaoptimset('MutationFcn', {@mutationuniform, mutationRate});
%
%   This will create an options structure specifying the mutation
%   function to be used is MUTATIONUNIFORM and the MUTATIONRATE is
%   user specified to be 0.05.
%

%   Copyright 2003-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/08/29 08:25:02 $

parents = obj.MutationParents;
mutationChildren = zeros(length(parents),obj.Dim);
for i=1:length(parents)
        child = obj.Population(parents(i),:);
        % Each element of the genome has mutationRate chance of being mutated.
        mutationPoints = find(rand(1,length(child)) < obj.MutationUniformRate);
        % each gene is replaced with a value chosen randomly from the range.
        range = obj.GlobalRange';
        % range can have one column or one for each gene.
        range = range(:,mutationPoints);
        lower = range(1,:);
        upper = range(2,:);
        span = upper - lower;
        child(mutationPoints) = lower + rand(1,length(mutationPoints)) .* span;
        mutationChildren(i,:) = child;
end

