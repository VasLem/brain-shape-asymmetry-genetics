classdef SM < mySuperClass
    % This is the abstract interface class for similarity measures
    % This class only interacts with transformation model objects
    % It is important for SM evaluations and Derivatives that the TM is
    % evaluated and Derivated beforhand in a main Routine;
    properties
        Target = [];% contains a handle to a target OBJECT;
        Correspondence = [];% Correspondences for Tmodel update, ICP based optimization, Correspondence is a meshObj or LMObj, similar to Tmodel.Evaluation;
    end
    properties (Abstract = true)
        % Attributes of these properties change in subclasses (ex. all
        % Dependent in combinedSM);
        Evaluation; % Tmodel point evaluations in similarity measure
        Derivative; % Tmodel point gradients in similarity measure        
        Fields;% The fields of LMobj or meshObj of which Similarity Measure is based on, mostly: Vertices
        nrSM;% mostly == 1, except for a combined SM
    end
    properties (Abstract = true, Dependent = true)
        TargetInfo;% DEPENDENT on target, retrieves the required information from the target;
    end
    properties (Dependent = true)
        nrFields;% number of Fields
        nrV;% number of points evaluated
    end
    methods %Constructor
        function obj = SM(varargin)
            obj = obj@mySuperClass(varargin{:});
        end
    end
    methods % Special Setting & Getting
        function nr = get.nrFields(obj)
            nr = length(obj.Fields);
        end
        function obj = set.Correspondence(obj,cor)
            if isempty(cor) 
               if ~isempty(obj.Correspondence), delete(obj.Correspondence); end 
               obj.Correspondence = [];
               return;
            end
            switch class(cor)
                case 'meshObj'
                    if isempty(obj.Correspondence), obj.Correspondence = fastClone(cor); return; end
                    if ~strcmp(class(obj.Correspondence),'meshObj'), delete(obj.Correspondence);obj.Correspondence = fastClone(cor); return; end
                    fastCopy(cor,obj.Correspondence);
                case 'LMObj'
                    if isempty(obj.Correspondence), obj.Correspondence = fastClone(cor); return; end
                    if ~strcmp(class(obj.Correspondence),'LMObj'), delete(obj.Correspondence);obj.Correspondence = fastClone(cor); return; end
                    fastCopy(cor,obj.Correspondence);
                otherwise
                   error('Correspondence Input is nor a meshObj, nor a LMobj'); 
            end
        end       
        function cor = get.Correspondence(obj)
                 cor = obj.Correspondence;
                 % Cannot call TM.validEval in this function, infinite Loop of recursions!!!
                 % isempty(obj.Correspondence) == ~validH(obj,'Correspondence');
                 if ~mySuperClass.isH(cor), cor = []; return; end
        end
        function target = get.Target(obj)
                 %disp('getting target');
                 target = obj.Target;
                 if ~mySuperClass.isH(target), target = []; end
        end
        function obj = set.Target(obj,target)
                 % allows you to make sure a correct target is given and
                 % that the target is properly preprocessed, e.g Target
                 % must always be a handles
                 %disp('setting target');
                 if ~mySuperClass.isH(target), return; end
                 target = preSetTarget(obj,target);
                 obj.Target = target; 
        end
        function nr = get.nrV(obj)
            nr = size(obj.Evaluation,2);
        end
    end
    methods % InterFace functions
        function delete(obj) %#ok<INUSD>
           if ~isempty(obj.Correspondence), delete(obj.Correspondence); end
           % Target is not deleted, Target is a kind of parent to SM
        end
        function struc = obj2struc(obj)
            % converting relevant information
            struc = obj2struc@mySuperClass(obj);
            struc.Evaluation = obj.Evaluation;
            struc.Derivative = obj.Derivative;
            struc.Correspondence = [];
            %struc.Target = [];
            if ~isempty(obj.Correspondence),struc.Correspondence = obj2struc(obj.Correspondence);end
            %if ~isempty(obj.Target), struc.Target = obj2struc(obj.Target);end
            struc.Fields = obj.Fields;
            struc.nrSM = obj.nrSM;
        end
        function obj = struc2obj(obj,struc)
            obj.Evaluation = struc.Evaluation;
            obj.Derivative = struc.Derivative;
            obj.Correpondence = [];
            if ~isempty(struc.Correspondence),obj.Correspondence = struc2obj(eval(struc.Correspondence.Type),struc.Correspondence);end
            %if ~isempty(struc.Target),obj.Target = struc2obj(eval(struc.Target.Type),struc.Target);end
            obj.Fields = struc.Fields;
        end
        function copy(obj,cobj)
                 copy@mySuperClass(obj,cobj);
                 cobj.Evaluation = obj.Evaluation;
                 cobj.Derivative = obj.Derivative;
                 cobj.Target = obj.Target;% Target is not cloned, Target is a kind of Parent of SM
                 if ~isempty(obj.Correspondence), cobj.Correspondence = fastClone(obj.Correspondence); end
                 %cobj.Fields = obj.Fields; Fields are initialized and are
                 %never changed
        end             
        function out = clear(obj)
                 if nargout == 0
                    obj = clone(obj);
                    out = obj;
                 end
                 obj.Evaluation = [];
                 obj.Derivative = [];
                 obj.Correspondence = [];
        end
        function out = initialize(obj,Tmodel)
                 if nargout == 1
                    obj = clone(obj);
                    out = obj;
                 end
                 if isempty(obj.Target), return; end
                 eval(obj,Tmodel);
        end
        function out = updateCorrespondence(obj,Tmodel) %#ok<INUSD>
                 warning([class(obj) ': is not able to update correspondences']); %#ok<WNTAG>
                 if nargout == 1, out = []; end
        end
        function out = isField(obj,field)
            if isempty(find(strcmp(obj.Fields,field),1)), out = false; return; end
            out = true;
        end
        function out = preSetTarget(obj,target) %#ok<INUSL>
            %disp('dummy preset target');
            out = target;%dummy
        end
        function out = sumEvaluation(obj,w)
            if nargin < 2
               out = sum(obj.Evaluation);
            elseif(numel(w)==1)
               out = w*sum(obj.Evaluation);
            else
               out = sum(w.*obj.Evaluation); 
            end
        end
        function out = sumDerivative(obj,w)
            if nargin < 2
               out = sum(obj.Derivative);
            elseif(numel(w)==1)
               out = w*sum(obj.Derivative);
            else
               out = sum(w.*obj.Derivative); 
            end
        end
    end
    methods (Abstract = true) %Abstract Interfaces Functions
            out = eval(obj,Tmodel);
            out = derivate(obj,Tmodel);
    end
end % classdef

