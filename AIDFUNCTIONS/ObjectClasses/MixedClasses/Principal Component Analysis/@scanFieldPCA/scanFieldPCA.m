classdef scanFieldPCA < scanPCA
    % This is an intermediate abstract class between scanPCA and every other
    % type of PCA class dealing with a single field of objects represented by meshes in particular faces
    % This class is storing methods that are applicable to the relevant
    % subclasses
    properties
        Average= [];% The average now is an object, having all the fields on which subclasses operate, the handling of the correct field is dealt with in the approriate subclasses;
    end
    methods % Constructor
        function obj = scanFieldPCA(varargin)
            obj = obj@scanPCA(varargin{:});         
        end
    end
    methods % Special Setting and Getting
        function out = get.Average(obj)
           out = obj.Average;
           if ~superClass.isH(out), out = [];return; end
        end
        function obj = set.Average(obj,in)
           if ~isempty(obj.Average)&&~(in==obj.Average), delete(obj.Average); end
           obj.Average = in; 
        end
    end
    methods
        function checkAverage(obj)
            if isempty(obj.Average),obj.Average = meshObj;end
        end
    end
end