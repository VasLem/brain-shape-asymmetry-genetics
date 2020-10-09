function mutationChildren = mutationgaussian(obj)
%MUTATIONGAUSSIAN Gaussian mutation.
%   MUTATIONCHILDREN = MUTATIONGAUSSIAN(PARENTS,OPTIONS,GENOMELENGTH,...
%   FITNESSFCN,STATE,THISSCORE,THISPOPULATION,SCALE,SHRINK) Creates the
%   mutated children using the Gaussian distribution.
%
%   SCALE controls what fraction of the gene's range is searched. A
%   value of 0 will result in no change, a SCALE of 1 will result in a
%   distribution whose standard deviation is equal to the range of this gene.
%   Intermediate values will produce ranges in between these extremes.
%
%   SHRINK controls how fast the SCALE is reduced as generations go by.
%   A SHRINK value of 0 will result in no shrinkage, yielding a constant search
%   size. A value of 1 will result in SCALE shrinking linearly to 0 as
%   GA progresses to the number of generations specified by the options
%   structure. (See 'Generations' in GAOPTIMSET for more details). Intermediate
%   values of SHRINK will produce shrinkage between these extremes.
%   Note: SHRINK may be outside the interval (0,1), but this is ill-advised.
%
%   Example:
%     options = gaoptimset('MutationFcn',{@mutationgaussian});
%
%   This specifies that the mutation function used will be
%   MUTATIONGAUSSIAN, and since no values for SCALE or SHRINK are specified
%   the default values are used.
%
%     scale = 0.5; shrink = 0.75;
%     options = gaoptimset('MutationFcn',{@mutationgaussian,scale,shrink});
%
%   This specifies that the mutation function used will be
%   MUTATIONGAUSSIAN, and the values for SCALE or SHRINK are specified
%   as 0.5 and 0.75 respectively.
%

%   Copyright 2003-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2011/05/09 00:50:35 $

% Use default parameters if the are not passed in.
% If these defaults are not what you prefer, you can pass in your own
% values when you set the mutation function:
%
% options.MutationFunction = { mutationgaussian, 0.3, 0} ;
%
parents = obj.Population(obj.MutationParents,:);
scales = obj.MutationScale*obj.FaceSigmas';
mutationChildren = parents + repmat(scales,size(parents,1),1).*randn(size(parents));
end
