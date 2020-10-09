classdef combinedSM < SM
    % This is the abstract interface class for combined similarity measures
    % multiple smeasures with regard to 1 target;
    % This class only interacts with transformation model objects
    % It is important for SM evaluations and Derivatives that the TM is
    % evaluated and Derivated beforhand in a main Routine;
    properties
        Smeasure = {};% List of similarity measures
        SmeasureW = [];% The weight (importance) of every Smeasure, default = 1 for every Smeasure (equally important)
    end
    properties (Dependent = true)
        Evaluation;% is now dependent on similarity measures part of the combined SM
        Derivative;% is now dependent on similarity measures part of the combined SM
        TargetInfo;% DEPENDENT on target, retrieves the required information from the target;
        Fields;% The fields of LMobj or meshObj of which Similarity Measure is based on, mostly: Vertices
        nrSM;% number of similarity measures
        sumSMW;% sum of similarity weights
    end
    methods %Constructor
        function obj = combinedSM(varargin)
            obj = obj@SM(varargin{:});
            if nargin > 0
             Input = find(strcmp(varargin, 'Smeasure'));if ~isempty(Input), obj.Smeasure = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'SmeasureW'));if ~isempty(Input), obj.SmeasureW = varargin{Input+1}; end
            end
        end
    end
    methods % Special Setting & Getting
        function nr = get.nrSM(obj)
            nr = length(obj.Smeasure);
        end
        function obj = set.nrSM(obj,nr) %#ok<INUSD>
            %dummy
        end
        function w = get.SmeasureW(obj)
                 w = obj.SmeasureW;
                 if isempty(w)
                    w = ones(1,obj.nrSM);% Defaults
                 end
        end
        function obj = set.SmeasureW(obj,w)
                 if numel(w) == 1
                    w = w*ones(1,obj.nrSM);
                 end
                 obj.SmeasureW = w;
        end
        function val = get.sumSMW(obj)
            val = sum(obj.SmeasureW);
        end
        function fields = get.Fields(obj)
                 if obj.nrSM==0, fields = {}; return; end
                 fields = {};
                 for i=1:1:obj.nrSM
                     fields = union(fields,obj.Smeasure{i}.Fields); %#ok<AGROW>
                 end            
        end
        function obj = set.Fields(obj,field) %#ok<INUSD>
                 % dummy
        end
        function val = get.Evaluation(obj)
                 if isempty(obj.Smeasure), val = []; return; end
                 val = zeros(obj.nrSM,size(obj.Smeasure{1}.Evaluation,2));
                 for i=1:1:obj.nrSM
                     if isempty(obj.Smeasure{i}.Evaluation), continue; end
                     val(i,:) = obj.SmeasureW(i)*obj.Smeasure{i}.Evaluation;
                 end
        end
        function obj = set.Evaluation(obj,val)
                  if isempty(obj.Smeasure), return; end
                  for i=1:1:obj.nrSM
                      obj.Smeasure{i}.Evaluation = val(i,:);
                  end
        end
        function val = get.Derivative(obj)
                 if isempty(obj.Smeasure), val = []; return; end
                 val = zeros(obj.nrSM,size(obj.Smeasure{1}.Derivative,2));
                 for i=1:1:obj.nrSM
                     if isempty(obj.Smeasure{i}.Derivative), continue; end
                     val(i,:) = obj.SmeasureW(i)*obj.Smeasure{i}.Derivative;
                 end
        end
        function obj = set.Derivative(obj,val)
                  if isempty(obj.Smeasure), return; end
                  for i=1:1:obj.nrSM
                      obj.Smeasure{i}.Derivative = val(i,:);
                  end
        end        
        function info = get.TargetInfo(obj)
            info = obj.Target;
        end
        function out = get.Smeasure(obj)
                 out = obj.Smeasure;
                 if isempty(out), return; end
                 index = (1:length(out));
                 badindex = [];
                 for i=1:1:length(out)
                     if ~mySuperClass.isH(out{i}), badindex = [badindex i]; end %#ok<AGROW>
                 end
                 restindex = setdiff(index,badindex);
                 out = out(restindex);
        end
        function obj = set.Smeasure(obj,in)
            if isempty(obj.Smeasure), obj.Smeasure = in; return; end
            for i=1:1:max(obj.nrSM,length(in))
                obj.Smeasure{i} = in{i};
            end
        end
    end
    methods % InterFace functions
        function delete(obj) %#ok<INUSD>
           delete@SM(obj)
           if isempty(obj.Smeasure), return; end
           tmp = obj.Smeasure;
            for i=1:1:obj.nrSM
                delete(tmp{i});
            end
        end
        function struc = obj2struc(obj)
            % converting relevant information
            struc = obj2struc@SM(obj);
            struc.SmeasureW = obj.SmeasureW;
            if isempty(obj.Smeasure), struc.Smeasure = {}; return; end
            for s=1:1:obj.nrSM
                struc.Smeasure{s} = obj2struc(obj.Smeasure{s});
            end
        end
        function obj = struc2obj(obj,struc)
            struc2obj@SM(obj,struc);
            obj.SmeasureW = struc.SmeasureW;
            if isempty(struc.Smeasure),obj.Smeasure = {}; return; end
            for s=1:1:length(struc.Smeasure)
                obj.Smeasure{s} = struc2obj(eval(struc.Smeasure{s}.Type),struc.Smeasure{s});
            end           
        end
        function copy(obj,cobj)
                 copy@SM(obj,cobj);
                 cobj.SmeasureW = obj.SmeasureW;
                 cobj.Smeasure = {};
                 if isempty(obj.Smeasure), return; end
                 for i=1:1:obj.nrSM
                     cobj.Smeasure{i} = clone(obj.Smeasure{i});
                 end
        end             
        function out = clear(obj)
                 if nargout == 0
                    obj = clone(obj);
                    out = obj;
                 end
                 clear@SM(obj);
                 if isempty(obj.Smeasure), return; end
                 for i=1:1:obj.nrSM
                     clear(obj.Smeasure{i});
                 end
        end
        function out = initialize(obj,Tmodel)
                 if nargout == 1
                    obj = clone(obj);
                    out = obj;
                 end
                 if isempty(obj.Target), return; end
                 eval(obj,Tmodel);
        end
        function out = preSetTarget(obj,target) %#ok<INUSL>
            out = target;%dummy
            if isempty(obj.Smeasure), return; end
            for i=1:1:obj.nrSM
                obj.Smeasure{i}.Target = target;
            end
        end
        function out = sumEvaluation(obj,w)
            if isempty(obj.Smeasure), out = []; return; end
            out = obj.Evaluation;
            if nargin < 2
               out = sum(sum(out));
            elseif(numel(w)==1)
               out = w*sum(sum(out));
            else
               out = sum(sum(repmat(w,obj.nrSM,1).*out));
            end
        end
        function out = sumDerivative(obj,w)
            if isempty(obj.Smeasure), out = []; return; end
            out = obj.Derivative;
            if nargin < 2
               out = sum(sum(out));
            elseif(numel(w)==1)
               out = w*sum(sum(out));
            else
               out = sum(sum(repmat(w,obj.nrSM,1).*out));
            end
        end
    end
end % classdef

