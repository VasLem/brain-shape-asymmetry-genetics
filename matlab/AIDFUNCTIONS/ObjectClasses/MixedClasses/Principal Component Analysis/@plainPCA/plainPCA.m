classdef plainPCA < PCA
    properties
        Average = [];% Average 
    end
    properties(Dependent = true)
        AvgVec;
    end
    methods % Constructor
        function obj = plainPCA(varargin)
            obj = obj@PCA(varargin{:});         
        end
    end
    methods % Special Setting & Getting
        function out = get.AvgVec(obj)
            out = obj.Average;
        end
        function obj = set.AvgVec(obj,in)
            obj.Average = in;
        end
    end
    methods % Specific Interface Functions
        function out = Struc2Vec(obj,in) %#ok<INUSL>
           % converts an Structure 2 a vector 
           out = in;
        end
        function out = Vec2Struc(obj,in) %#ok<INUSL>
            % converts a vector 2 structure
            out = in;
        end
        function out = IndexStruc2Vec(obj,in) %#ok<INUSL>
            % converts an index defined on a structure to an index defined
            % on a vector;
            out = in;
        end
        function out = IndexVec2Struc(obj,in) %#ok<INUSL>
            % converts an index defined on a vector to an index defined on
            % a structure
            out = in;
        end
        function out = WeightStruc2Vec(obj,in) %#ok<INUSL>
            % converts weights defined on a structure into weights defined
            % on a vector
            out = in;
        end
        function out = WeightVec2Struc(obj,in) %#ok<INUSL>
            % converts weights defined on a vector into weights defined on
            % a structure;
            out = in;
        end
        function out = getData(obj,in) %#ok<INUSL>
           out = in; 
        end
    end
end