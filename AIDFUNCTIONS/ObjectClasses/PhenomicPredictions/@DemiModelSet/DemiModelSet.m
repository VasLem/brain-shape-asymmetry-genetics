classdef DemiModelSet < superClass
    % This class will is a CONTAINER class to accumulate multiple demiModels into a set, and
    % handles the input parsing, accumulated error scoring...
    % within the EvoMorph algorithm this provides the fitness function
    properties
        DM;% list of demimodels
        Input;% of a face to reconstruct
        ID;% ID of the face to reconstruct
        DMWeights = [];% individual weights for demimodels
    end
    properties (Dependent = true)
        nrModels;% number of demimodels
        toUse;% the demimodels to use (is dependent on the avalaible information in input values)
        nrUse;% nr of demimodels to use
        DMTags;% list of DM tags
        DMLevels;% list of DM Levels
    end
    methods % Constructor
        function obj = DemiModelSet(varargin)
            obj = obj@superClass(varargin{:});         
        end
    end
    methods % Special Setting & Getting
        function out = get.nrModels(obj)
            out = length(obj.DM);
        end
        function out = get.toUse(obj)
            out = find(~isnan(obj.Input));
        end
        function out = get.nrUse(obj)
            out = length(obj.toUse);
        end
        function out = get.DMWeights(obj)
            out = obj.DMWeights;
            if isempty(out), out = ones(1,obj.nrModels); end
        end
        function out = get.DMTags(obj)
           out = cell(1,obj.nrModels);
           for i=1:1:obj.nrModels
               out{i} = obj.DM{i}.Tag;
           end
        end
        function out = get.DMLevels(obj)
           out = cell(1,obj.nrModels);
           for i=1:1:obj.nrModels
               out{i} = obj.DM{i}.Level;
           end
        end
    end
    methods % General Interface Functions        
        function out = getScore(obj,in)
            n = size(in,2);
            scores = zeros(n,obj.nrUse);
            parfor i=1:obj.nrUse
                scores(:,i) = getScore(obj.DM{obj.toUse(i)},in,obj.Input(obj.toUse(i)));%#ok<*PFBNS>toc;
            end
            out = sum(scores.*repmat(obj.DMWeights(obj.toUse),n,1),2)';%./sum(obj.DMWeights(obj.toUse));% this is a straigthforward weighted summing of errors 
        end
        function parseInput(obj,in)
           [~,~,raw] = xlsread(in);
           ind = find(strcmp('ID',raw(:,1)));
           if ~isempty(ind)
               obj.ID = raw{ind,2};
           else
               disp('ID is missing');
           end
           val = nan*zeros(1,obj.nrModels);
           for i=1:obj.nrModels
              ind =  find(strcmp(obj.DM{i}.Tag,raw(:,1)));
              if~isempty(ind) 
                val(i) = raw{ind,2};
              end
           end
           obj.Input = val;
        end
        function relateWeights2Acc(obj)
           % this is especially needed for genotype classificatin
           out = ones(1,obj.nrModels);
           for i=1:1:obj.nrModels
               switch obj.DM{i}.XType
                   case 'Categorical'
                       out(i) = obj.DM{i}.Acc/100;
                   case 'Continous'
                       continue;
               end
           end
           obj.DMWeights = out;
        end
    end
end