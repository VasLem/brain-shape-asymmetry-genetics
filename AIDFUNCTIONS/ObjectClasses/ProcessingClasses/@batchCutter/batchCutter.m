classdef batchCutter < batchConvertor
    properties
       Index = [];
    end
    methods %Constructor
        function obj = batchCutter(varargin)
          obj = obj@batchConvertor(varargin{:});
          obj.InputFormat = 'Mat Object';
          obj.OutputFormat = 'Mat Object';
          obj.Overwrite = false;
        end
    end
    methods% Interface Functions
        function scan = innerFunction(obj,scan)
            crop(scan,'VertexIndex',obj.Index);
        end     
    end
end