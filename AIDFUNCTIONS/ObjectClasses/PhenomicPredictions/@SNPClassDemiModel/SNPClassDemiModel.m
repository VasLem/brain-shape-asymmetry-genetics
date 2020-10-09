classdef SNPClassDemiModel < superClass
    % This is a specific class for SNP classification, and uses multiple
    % TwoClassDemiModel objects
   properties
       TexLab = {'AA' 'AB' 'BB'};
       CatLab = [-1 0 1];
       AAABModel = TwoClassDemiModel;
       AABBModel = TwoClassDemiModel;
       ABBBModel = TwoClassDemiModel;
   end
   methods% constructor
        function obj = SNPClassDemiModel(varargin)
            obj = obj@superClass(varargin{:});         
        end
   end
   methods % special setting and getting
   end
   methods % general interface functions
   end
end