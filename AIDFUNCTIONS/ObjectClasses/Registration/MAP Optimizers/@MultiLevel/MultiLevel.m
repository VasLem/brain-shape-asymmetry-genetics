classdef MultiLevel < mySuperClass
    % This is the abstract interface class for Multi Level (Hierachical)
    % Shemes, affecting the used number of Points in a MAP objective
    % function optimization
    properties
        MinPoints = 1000;% minimum amounts of Points used, lowest level
        MaxPoints = [];% Maximum amount of Points used, highest level,if [] == floating.nrV
        nrLevels = 4;% Number of Levels in the Sheme
        Level = 1;% Current Level
    end
    properties (Dependent = true)
        nrPoints;
        PointFraction;
    end
    methods %Constructor
        function obj = MultiLevel(varargin)
          obj = obj@mySuperClass(varargin{:});
          if nargin > 0
             Input = find(strcmp(varargin, 'MinPoints'));if ~isempty(Input), obj.MinPoints = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'MaxPoints'));if ~isempty(Input), obj.MaxPoints = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'nrLevels'));if ~isempty(Input), obj.nrLevels = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'Level'));if ~isempty(Input), obj.Level = varargin{Input+1}; end
          end
        end
    end
   methods % Special Getting & Setting
       function out = get.PointFraction(obj)
           out = (obj.MaxPoints-obj.MinPoints)/(obj.nrLevels-1); 
       end
       function out = get.nrPoints(obj)
           out = ceil(obj.MinPoints + (obj.Level-1)*obj.PointFraction);
       end
   end
   methods % InterFace functions
        function copy(obj,cobj)
                 copy@mySuperClass(obj,cobj);
                 cobj.MinPoints = obj.MinPoints;
                 cobj.MaxPoints = obj.maxPoints;
                 cobj.NrLevels = obj.nrLevels;
                 cobj.Level = obj.Level;
        end
        function struc = obj2struc(obj)
                 struc = obj2struc@mySuperClass(obj);
                 struc.MinPoints = obj.MinPoints;
                 struc.MaxPoints = obj.maxPoints;
                 struc.NrLevels = obj.nrLevels;
                 struc.Level = obj.Level;
        end
        function obj = struc2obj(obj,struc)
                 obj.MinPoints = struc.MinPoints;
                 obj.MaxPoints = struc.maxPoints;
                 obj.NrLevels = struc.nrLevels;
                 obj.Level = struc.Level;
        end
        function str = strInfo(obj)
            str = ['ML(' num2str(obj.Level) ',' num2str(obj.nrLevels) ',' num2str(obj.nrPoints) ')/'];
        end
        function displayInfo(obj)
            display(strInfo(obj));
        end
   end
end % classdef