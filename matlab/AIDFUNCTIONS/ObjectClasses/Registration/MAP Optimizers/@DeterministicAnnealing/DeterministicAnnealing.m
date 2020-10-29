classdef DeterministicAnnealing < superClass
    % This is the abstract interface class for Deterministic Anealing
    % Shemes, affecting the Temperature parameter in a MAP objective
    % function
    properties
        StartTemp = 1;% initial Temperature
        TempFraction = 500;% TempFraction of the initial Temperature
        Rate = 0.9;% Annealing rate
        Step = 1; % DeterministicAnnealing iteration (step)
        MaxIter = 3;% Innerloop iterations are typically less, when using DeterministicAnnealing
        FinalZero = true;
        %FinalTemp = [];
    end
    properties (Dependent = true)
        FinalTemp; % Final Temperature
        nrSteps;% Number of annealing steps
    end
    methods %Constructor
        function obj = DeterministicAnnealing(varargin)
          obj = obj@superClass(varargin{:});
        end
    end
   methods % Special Getting & Setting
       function out = get.FinalTemp(obj)
           %if ~isempty(obj.FinalTemp), out = obj.FinalTemp; return; end
           %if obj.FinalZero, out = 0; return; end
           out = obj.StartTemp/obj.TempFraction;
       end
       function obj = set.FinalTemp(obj,in)
           %obj.FinalTemp = in;
           obj.TempFraction = obj.InitTemp/in;
       end
       function out = get.nrSteps(obj)
           out = ceil(log(obj.FinalTemp/obj.StartTemp)/log(obj.Rate));
       end
       function obj = set.nrSteps(obj,in)
           t = obj.StartTemp;
           for i=1:1:in
               t = t*obj.Rate;
           end
           obj.FinalTemp = t;
       end
   end
   methods % InterFace functions
        function str = strInfo(obj)
            str = ['DA(' num2str(obj.Step) ',' num2str(obj.nrSteps) ')/'];
        end
        function displayInfo(obj)
            display(strInfo(obj));
        end
        function updateTemp(obj,MapFunction)
            if obj.FinalZero&&obj.Step==obj.nrSteps, MapFunction.Temp = 0; return; end
            MapFunction.Temp = obj.Rate*MapFunction.Temp;         
        end
   end
end % classdef