classdef singleCP < CP
    % This is the abstract interface class Complete Processes
    properties
        InlierP = []; % link to Inlier Process
        OutlierP = []; % link to Outlier Process
        LatentV = []; % link to Latent variable
    end
    properties (Dependent = true)
        Target;% Target handle is stored in CompleteProcess.InlierP.Smeasure.Target;
        Smeasure;% Smeasure handle is stored in CompleteProcess.InlierP
        Kappa;% Statistical significance of Inlierp vs Outlierp
        Fields;
        B;% B values = LatentV.Value
        Correspondence;% Dependent on Correspondence stored in CompleteProcess.InlierProces.Smeasure.Correspondence
        Temp;
        FastEval;% Fast (re)-evaluation of inlier process
    end
    methods %Constructor
        function obj = singleCP(varargin)
          obj = obj@CP(varargin{:});
        end
    end  
    methods % Special Setting and Getting
        function obj = set.InlierP(obj,input)
                 if~isempty(obj.InlierP), delete(obj.InlierP); end
                 if isempty(input), obj.InlierP = []; return; end
                 obj.InlierP = input;
                 input.CompleteP = obj;
        end
        function out = get.InlierP(obj)
                 out = obj.InlierP;
                 if ~mySuperClass.isH(out), out = []; end
        end
        function obj = set.OutlierP(obj,input)
                 if~isempty(obj.OutlierP), delete(obj.OutlierP); end
                 if isempty(input), obj.OutlierP = []; return; end
                 obj.OutlierP = input;
                 input.CompleteP = obj;
        end
        function out = get.OutlierP(obj)
                 out = obj.OutlierP;
                 if ~mySuperClass.isH(out), out = []; end
        end
        function obj = set.LatentV(obj,input)
                 if~isempty(obj.LatentV), delete(obj.LatentV); end
                 if isempty(input), obj.LatentV = []; return; end
                 obj.LatentV = input;
                 input.CompleteP = obj;
        end
        function out = get.LatentV(obj)
                 out = obj.LatentV;
                 if ~mySuperClass.isH(out), out = []; end
        end
        function obj = set.Target(obj,target)
                 if isempty(obj.InlierP), return; end
                 obj.InlierP.Target = target;
        end
        function target = get.Target(obj)
            if isempty(obj.InlierP), target = []; return; end
            target = obj.InlierP.Target;
        end
        function obj = set.Smeasure(obj,smeasure)
                 if isempty(obj.InlierP), return; end
                 obj.InlierP.Smeasure = smeasure;
        end
        function smeasure = get.Smeasure(obj)
            if isempty(obj.InlierP), smeasure = []; return; end
            smeasure = obj.InlierP.Smeasure;
        end
        function obj = set.Kappa(obj,in)
                 if isempty(obj.InlierP), return; end
                 obj.InlierP.Kappa = in;
        end
        function out = get.Kappa(obj)
            if isempty(obj.InlierP), out = []; return; end
            out = obj.InlierP.Kappa;
        end
        function obj = set.Correspondence(obj,cor)
                 if isempty(obj.InlierP), return; end
                 obj.InlierP.Correspondence = cor;
                 %InlierP Already has a link to correspondence stored in Smeasure
        end
        function cor = get.Correspondence(obj)
            if isempty(obj.InlierP), cor = []; return; end
            cor = obj.InlierP.Correspondence;
        end
        function B = get.B(obj)
            if isempty(obj.LatentV), B = []; return; end
            B = obj.LatentV.Value;
        end
        function fields = get.Fields(obj)
            if isempty(obj.InlierP), fields = {}; return; end
            fields = obj.InlierP.Fields;
        end
        function obj = set.Fields(obj,field) %#ok<INUSD>
                 % dummy
        end
        function obj = set.Temp(obj,in)
                 if isempty(obj.InlierP), return; end
                 obj.InlierP.Temp = in;
        end
        function out = get.Temp(obj)
            if isempty(obj.InlierP), out = []; return; end
            out = obj.InlierP.Temp;
        end
        function obj = set.FastEval(obj,in)
                 if isempty(obj.InlierP), return; end
                 obj.InlierP.FastEval = in;
        end
        function out = get.FastEval(obj)
            if isempty(obj.InlierP), out = []; return; end
            out = obj.InlierP.FastEval;
        end
    end       
    methods % InterFace functions 
        function copy(obj,cobj)
                 copy@CP(obj,cobj);
                 if ~isempty(obj.InlierP),cobj.InlierP = clone(obj.InlierP);end
                 if ~isempty(obj.OutlierP),cobj.OutlierP = clone(obj.OutlierP);end
                 if ~isempty(obj.LatentV),cobj.LatentV = clone(obj.LatentV);end          
        end
        function delete(obj) 
            if ~isempty(obj.InlierP),delete(obj.InlierP);end
            if ~isempty(obj.OutlierP),delete(obj.OutlierP);end
            if ~isempty(obj.LatentV),delete(obj.LatentV);end
        end
        function out = clear(obj)
                 if nargout == 1
                    obj = clone(obj);
                    out = obj;
                 end
                 clear@CP(obj);
                 if ~isempty(obj.InlierP),clear(obj.InlierP);end
                 if ~isempty(obj.OutlierP),clear(obj.OutlierP);end
                 if ~isempty(obj.LatentV),clear(obj.LatentV);end
        end
        function struc = obj2struc(obj)
            % converting relevant information
            struc = obj2struc@CP(obj);
            struc.InlierP = [];struc.OutlierP = [];struc.LatentV = [];
            if ~isempty(obj.InlierP),  struc.InlierP = obj2struc(obj.InlierP);end
            if ~isempty(obj.OutlierP), struc.OutlierP = obj2struc(obj.OutlierP);end
            if ~isempty(obj.LatentV),  struc.LatentV = obj2struc(obj.LatentV);end
        end
        function obj = struc2obj(obj,struc)
                 obj = struc2obj@CP(obj,struc);
                 obj.InlierP = [];obj.OutlierP = [];obj.LatentV = [];
                 if ~isempty(struc.InlierP), obj.InlierP = struc2obj(eval(struc.InlierP.Type),struc.InlierP); end
                 if ~isempty(struc.OutlierP), obj.OutlierP = struc2obj(eval(struc.OutlierP.Type),struc.OutlierP); end
                 if ~isempty(struc.LatentV), obj.LatentV = struc2obj(eval(struc.LatentV.Type),struc.LatentV); end 
        end
        function out = eval(obj,Tmodel)
                 %if isempty(obj.InlierP), error('Cannot eval CP: InlierP is not set'); end
                 eval(obj.InlierP,Tmodel);
                 eval(obj.OutlierP,Tmodel);
                 tmpout = obj.InlierP.Evaluation + obj.OutlierP.Evaluation;
                 if nargout == 1, out = tmpout; return; end
                 obj.Evaluation = tmpout;          
        end
        function out = derivate(obj,Tmodel)
                 derivate(obj.InlierP,Tmodel);
                 derivate(obj.OutlierP,Tmodel);
                 tmpout = obj.InlierP.Derivative + obj.OutlierP.Derivative;
                 if nargout == 1, out = tmpout; return; end
                 obj.Derivative = tmpout;
        end
        function out = update(obj,Tmodel)
                 if nargout == 1
                    obj = clone(obj);
                    out = obj;
                 end
                 update(obj.InlierP,Tmodel);% First Inlier Process update (similarity measures and parameters)
                 update(obj.OutlierP,Tmodel);% Second Outlier Process update (parameters dependent on Inlier process)
                 update(obj.LatentV,Tmodel);% Finnally reestimation of latent variables given new inlier and outlier processes
                 obj.FastEval = true;eval(obj,Tmodel);obj.FastEval = false;% Fast (re)-evaluation using new latentV values
        end
        function out = initialize(obj,Tmodel)
            if nargout == 1
               obj = clone(obj);
               out = obj;
            end
            initialize(obj.LatentV,Tmodel);% First setting LatentVariables equal to 1
            initialize(obj.InlierP,Tmodel);% Then update (initialize) Inlier Process
            initialize(obj.OutlierP,Tmodel);% Then update (initialize) Outlier Process
            obj.FastEval = true;eval(obj,Tmodel);obj.FastEval = false;% First (fast) evaluations
        end
        function out = updateCorrespondence(obj,Tmodel)
            % updating correspondences for ICP based optimization
            if isempty(obj.InlierP), error('Cannot update correspondences CP, no InlierP set'); end
            if nargout == 1
                out = updateCorrespondence(obj.InlierP,Tmodel);
            else
                updateCorrespondence(obj.InlierP,Tmodel);
            end           
        end
    end
end % classdef

