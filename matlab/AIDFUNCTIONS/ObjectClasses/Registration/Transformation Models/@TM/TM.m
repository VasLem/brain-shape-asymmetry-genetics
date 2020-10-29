classdef TM < superClass%TM < mySuperClass
    % This is the abstract interface class for Transformation models
    properties
        Evaluation = [];% floating surface points Evaluation in transformation model; This field is an OBJECT!!!
        Derivative = [];% floating surface points Derivative towards Active parameter in transformation model; THis property is a structure!!!!
        LnormEvaluation = [];% evaluation of the (log) Lnorm using current parameter setting
        LnormDerivative = [];% Derivative of the (log) Lnorm Gradient towards Active parameter
        Nu = 1;% user defined Lnorm weigthing factor (acts like regularization)
        ActiveField = [];% The Field of LMObj or meshObj on which the TM is currently active
        ActiveP = 1;% Active Parameter
        Temp = 0;% Temperature Parameter, regularizing Transformation model if ~==0 (0<Temp<1)
    end
    properties (Abstract = true)
        P;% Transformation parameters
        d;% Parameter deviations to compute numerical gradients
        Fields; % The fields of LMobj or meshObj on which the transformation applies, mostly: Vertices
    end
    properties (Dependent = true)
        nrP;% number of parameters
        nrV;% number of evaluated points
        nrFields;% number of Fields
        VertDerivative;% Point location(Vertices) derivative (stored in Derivative)
        TextDerivative;% Point Texture derivative (stored in Derivative)
        DistDerivative;% Point distance derivative (stored in Derivative)
    end
    methods %Constructor
        function obj = TM(varargin)
          obj = obj@superClass(varargin{:});

        end
    end  
    methods % Special Setting and Getting
        function obj = set.Derivative(obj,grad)
                 if (isempty(grad)&&isempty(obj.ActiveField)), obj.Derivative = []; return; end
                 %if isempty(obj.ActiveField), return; end
                 if ~isField(obj,obj.ActiveField), return; end%warning('No proper Active Field: action ignored'); return; end %#ok<WNTAG>
                 warning('off','All');
                 obj.Derivative.(obj.ActiveField) = grad;
                 warning('on','All');
        end
        function val = get.Derivative(obj)
                 if isempty(obj.ActiveField), val = obj.Derivative; return; end
                 if isfield(obj.Derivative,obj.ActiveField)
                    val = obj.Derivative.(obj.ActiveField);
                 else
                    val = [];
                 end
        end
        function grad = get.VertDerivative(obj)
                 old = obj.ActiveField;
                 obj.ActiveField = 'Vertices';
                 grad = obj.Derivative;
                 obj.ActiveField = old;
        end
        function obj = set.VertDerivative(obj,grad)
                 old = obj.ActiveField;
                 obj.ActiveField = 'Vertices';
                 obj.Derivative = grad;
                 obj.ActiveField = old;
        end
        function grad = get.TextDerivative(obj)
                 old = obj.ActiveField;
                 obj.ActiveField = 'TextureColor';
                 grad = obj.Derivative;
                 obj.ActiveField = old;
        end
        function obj = set.TextDerivative(obj,grad)
                 old = obj.ActiveField;
                 obj.ActiveField = 'TextureColor';
                 obj.Derivative = grad;
                 obj.ActiveField = old;
        end
        function grad = get.DistDerivative(obj)
                 old = obj.ActiveField;
                 obj.ActiveField = 'Distance';
                 grad = obj.Derivative;
                 obj.ActiveField = old;
        end
        function obj = set.DistDerivative(obj,grad)
                 old = obj.ActiveField;
                 obj.ActiveField = 'Distance';
                 obj.Derivative = grad;
                 obj.ActiveField = old;
        end
        function nr = get.nrP(obj)
            nr = length(obj.P);
        end
        function nr = get.nrV(obj)
            if isempty(obj.Evaluation), nr = 0; return; end
            nr = obj.Evaluation.nrV;
        end
        function nr = get.nrFields(obj)
            nr = length(obj.Fields);
        end
        function obj = set.Evaluation(obj,fs)
            if isempty(fs) 
               if ~isempty(obj.Evaluation), delete(obj.Evaluation); end 
               obj.Evaluation = [];
               return;
            end
            switch class(fs)
                case 'meshObj'
                    if isempty(obj.Evaluation), obj.Evaluation = fastClone(fs); return; end
                    if ~strcmp(class(obj.Evaluation),'meshObj'), delete(obj.Evaluation);obj.Evaluation = fastClone(fs); return; end
                    fastCopy(fs,obj.Evaluation);
                case 'LMObj'
                    if isempty(obj.Evaluation), obj.Evaluation = fastClone(fs); return; end
                    if ~strcmp(class(obj.Evaluation),'LMObj'), delete(obj.Evaluation);obj.Evaluation = fastClone(fs); return; end
                    fastCopy(fs,obj.Evaluation);
                otherwise
                   error('Evaluation Input is nor a meshObj, nor a LMobj'); 
            end
        end
        function eval = get.Evaluation(obj)
                 %if~isempty(obj.Evaluation), eval = []; return; end
                 % Cannot call mySuperClass.validH in this function,
                 % infinite Loop of recursions!!!
                 eval = obj.Evaluation;
                 % isempty(obj.Evaluation) == ~validH(obj,'Evaluation') or
                 % ~isH(obj.Evaluation)
                 if ~mySuperClass.isH(eval), eval = []; end                 
        end
        function obj = set.ActiveP(obj,index)
            if (index<1)||(index>obj.nrP), error('Active Parameter is not within range of TM parameters'); end
            obj.ActiveP = index;
        end
        function obj = set.ActiveField(obj,field)
                 obj.ActiveField = field;
                 postActiveField(obj,field);
        end
        function obj = set.Temp(obj,val)
            if (val<0)||(val>1), error('Temperature must be between 0 and 1'); end
            obj.Temp = val;
            postTemp(obj,val);
        end
    end       
    methods % InterFace functions
        function out = isField(obj,field)
            if isempty(find(strcmp(obj.Fields,field),1)), out = false; return; end
            out = true;
        end
        function diff = difference(obj,obj2)
                 if ~strcmp(obj.Type,obj2.Type), error('TMs are different, cannot take difference'); end
                 diff = norm(obj.P-obj2.P);
        end
        function out = transform(obj,p)
            switch class(p)
                case {'meshObj' 'LMObj' 'floatingS'}
                    eval(obj,p);
                    if nargout == 1
                       p = clone(p);
                       out = p;
                    end
                    fastCopy(obj.Evaluation,p);
                    updateChildren(p,'Transform',obj);
                    %delete(tmpout);
                otherwise
                    if ~size(p,1)==3, p = p';end
                    eval(obj,p);% out becomes a LMObj
                    out = obj.Evaluation.Vertices;
                    %delete(tmpout);
            end
        end
        function out = initialize(obj,p)
            if nargout == 1,obj = clone(obj);out = obj;end
            %p = TM.convertInput(p);
            obj.ActiveField = [];
            obj.P = zeros(obj.nrP,1);
            obj.Nu = 1;
            obj.Evaluation = p;% performs only a fastCopy
            copy(p,obj.Evaluation);% hence followed by a full copy (rendering becomes the same)
            obj.Derivative = [];
            obj.LnormEvaluation = 0;% cost for a transformation with P=zeros == 0;
            obj.LnormDerivative = [];    
        end
        function out = refresh(obj,fs)
            if nargout == 1,obj = clone(obj);out = obj;end
            clear(obj);
            eval(obj,fs);
            derivate(obj,fs);
            logLnorm(obj);
            logLnormDerivative(obj);
        end
        function out = clear(obj)
                 % Clearing the eval and deriv fields... 
                 if nargout == 1,obj = clone(obj);out = obj;end
                 obj.ActiveField = [];
                 obj.Evaluation = [];
                 obj.Derivative = [];
                 obj.LnormEvaluation = [];
                 obj.LnormDerivative = [];
        end
        function delete(obj) %#ok<INUSD>
           if ~isempty(obj.Evaluation), delete(obj.Evaluation); end
        end
        function out = chainRule(obj,in,p)%#ok<INUSD,INUSL>
            % needed to riple back pointGradients in combined transformtions
            % This is a dummy chain rule (out = in), if required needs to be
            % reimplemented in relevant subclasses
            out = in;
        end
        function postActiveField(obj,field) %#ok<INUSD>
                 return;
                 % Dummy function, only relevant in combinedTM
        end
        function postTemp(obj,val) %#ok<INUSD>
                 return;
                 % Dummy function, only relevant in combinedTM
        end
        function varargout = derivateField(obj,p,field)
            if iscell(field)
                nrF = length(field);
                if nrF==1, field = field{1}; end
            else
                nrF = 1;
            end
            old = obj.ActiveField;
            if nrF==1
                 switch field
                     case 'All'% Evaluate all fields given in Tmodel
                         for i=1:1:obj.nrFields
                             obj.ActiveField = obj.Fields{i};
                             if nargout >= 1
                                if i>nargout, return;end
                                varargout{i} = derivate(obj,p);
                             else
                                derivate(obj,p);
                             end
                         end
                     otherwise
                          if~isField(obj,field)
                             if nargout>=1, varargout{1} = []; end
                             return; 
                          end
                          obj.ActiveField = field;
                          if nargout>=1
                             varargout{1} = derivate(obj,p);
                          else
                             derivate(obj,p);
                          end
                 end
            elseif nrF>1% typically given by fields defined in SM/InlierProcess
                for i=1:1:nrF
                    obj.ActiveField = field{i};
                    if nargout >= 1
                       if i>nargout, return;end
                       varargout{i} = derivate(obj,p);
                    else
                       derivate(obj,p);
                    end
                end
            else
               error('Wrong Field input');
            end
            obj.ActiveField = old;
        end
        function updateLevel(obj,p,q,type) %#ok<INUSD>
            %dummy
        end
    end
    methods (Static = true) % Static interface functions
        function input = convertInput(input)
            switch class(input)
                case {'meshObj' 'LMObj'}
                    input = fastClone(input);
                case {'double' 'single'}
                    if ~size(input,1) == 3, input = input';end
                    if ~size(input,1) == 3; error('Input must be meshObj, LMObj or a list of 3D points'); end
                    input = LMObj('Vertices',input);
                otherwise
                    error('Input must be meshObj, LMObj or a list of 3D points');
            end
        end
        function p = getPoints(p)
           % function to get the point list from objects or pass on given point list
             switch class(p)
                 case {'meshObj' 'LMObj' 'floatingS'} 
                   p = p.Vertices;
                 otherwise
                  if ~size(p,1)==3, p = p';end
            end
        end
        function nr = getNrPoints(p)
           % function to get the point list from objects or pass on given point list
             p = TM.getPoints(p);
             nr = size(p,2);
        end
    end
    methods (Abstract = true) %Abstract Interfaces Functions
            out = eval(obj,p);
            out = match(obj,p,q,w);
            out = derivate(obj,p);% field: Vertices, TextureColor, Distance,...
            out = logLnorm(obj);
            out = logLnormDerivate(obj);
    end
end % classdef

%         function copy(obj,cobj)
%                  copy@mySuperClass(obj,cobj);
%                  cobj.P = obj.P;
%                  cobj.d = obj.d;
%                  if ~isempty(obj.Evaluation)
%                     cobj.Evaluation = clone(obj.Evaluation);
%                  else
%                     cobj.Evaluation = [];
%                  end
%                  cobj.VertDerivative = obj.VertDerivative;
%                  cobj.TextDerivative = obj.TextDerivative;
%                  cobj.DistDerivative = obj.DistDerivative;
%                  cobj.LnormEvaluation = obj.LnormEvaluation;
%                  cobj.LnormDerivative = obj.LnormDerivative;
%                  cobj.Nu = obj.Nu;
%                  cobj.Temp = obj.Temp;
%         end

%           if nargin > 0
%              Input = find(strcmp(varargin, 'P'));if ~isempty(Input), obj.P = varargin{Input+1}; end
%              Input = find(strcmp(varargin, 'd'));if ~isempty(Input), obj.d = varargin{Input+1}; end
%              Input = find(strcmp(varargin, 'Evaluation'));if ~isempty(Input), obj.Evaluation = varargin{Input+1}; end
%              Input = find(strcmp(varargin, 'Derivative'));if ~isempty(Input), obj.Derivative = varargin{Input+1}; end
%              Input = find(strcmp(varargin, 'LnormEvaluation'));if ~isempty(Input), obj.LnormEvaluation = varargin{Input+1}; end
%              Input = find(strcmp(varargin, 'LnormDerivative'));if ~isempty(Input), obj.LnormDerivative = varargin{Input+1}; end
%              Input = find(strcmp(varargin, 'Nu'));if ~isempty(Input), obj.Nu = varargin{Input+1}; end
%              Input = find(strcmp(varargin, 'Temp'));if ~isempty(Input), obj.Temp = varargin{Input+1}; end
%           end

%         function struc = obj2struc(obj)
%             % converting relevant information
%             struc = obj2struc@mySuperClass(obj);
%             struc.P = obj.P;
%             struc.d = obj.d;
%             if ~isempty(obj.Evaluation)
%                struc.Evaluation = obj2struc(obj.Evaluation);
%             else
%                struc.Evaluation = [];
%             end
%             struc.VertDerivative = obj.VertDerivative;
%             struc.TextDerivative = obj.TextDerivative;
%             struc.DistDerivative = obj.DistDerivative;
%             struc.LnormEvaluation = obj.LnormEvaluation;
%             struc.LnormDerivative = obj.LnormDerivative;
%             struc.Nu = obj.Nu;
%             struc.Temp = obj.Temp;
%         end
%         function obj = struc2obj(obj,struc)
%                  obj.P = struc.P;
%                  obj.d = struc.d;
%                  if ~isempty(struc.Evaluation)
%                     obj.Evaluation = struc2obj(eval(struc.Evaluation.Type),struc.Evaluation);
%                  else
%                     obj.Evaluation = [];
%                  end
%                  obj.VertDerivative = struc.VertDerivative;
%                  obj.TextDerivative = struc.TextDerivative;
%                  obj.DistDerivative = struc.DistDerivative;
%                  obj.LnormEvaluation = struc.LnormEvaluation;
%                  obj.LnormDerivative = struc.LnormDerivative;
%                  obj.Nu = struc.Nu;
%                  obj.Temp = struc.Temp;
%         end

