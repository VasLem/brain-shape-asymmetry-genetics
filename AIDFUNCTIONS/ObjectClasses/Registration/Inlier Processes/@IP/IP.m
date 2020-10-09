classdef IP < mySuperClass
    % This is the abstract interface class for Inlier Processes
    % These kind of processes, encapsulate similarity measures in a pdf,
    % increasing noise robustness and combined with an Outlier process into a completeProcess,
    % also increasing outlier robustness
    properties 
        CompleteP = []; % link to Parent, being the complete Process
        Smeasure = []; % similarity measures contained in Inlier Process
        Evaluation = [];% Inlier Process negative log evaluations over all floating surface points
        Derivative = [];% Inlier Process gradients in floating surface points
        PdfEvaluation = [];% Inlier Process pdf evaluation of floating surface points
        Kappa = 3;% Statistical significance level, being regarded to be inliers
        Temp = 0;% Temperature parameter, relaxes the InlierProcess, Affects the Parameter estimation (0<Temp<1)
        FastEval = false;% If set to true, evaluation is performed without re-evaluating similarity measures
    end
    properties (Abstract = true)
        PdfP;% Inlier Process parameters
    end
    properties (Dependent = true)
        nrPdfP;% number of parameters
        nrSM;% number of similarity measures
        nrV;% number of points evaluated
        Fields% Dependent on similarity measures, Fields on which the IP operates
        Target;% Dependent on Similarity measures.
        nrFields;% number of Fields
        Correspondence = [];% Correspondences for Tmodel update, ICP based optimization, Correspondence is a meshObj or LMObj, similar to Tmodel.Evaluation;and is combination of Correspondences obtained in the similarity measures
    end
    methods %Constructor
        function obj = IP(varargin)
          obj = obj@mySuperClass(varargin{:});
          if nargin > 0
             Input = find(strcmp(varargin, 'PdfP'));if ~isempty(Input), obj.PdfP = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'CompleteP'));if ~isempty(Input), obj.CompleteP = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'Smeasure'));if ~isempty(Input), obj.Smeasure = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'Evaluation'));if ~isempty(Input), obj.Evaluation = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'Derivative'));if ~isempty(Input), obj.Derivative = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'PdfEvaluation'));if ~isempty(Input), obj.PdfEvaluation = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'Correspondence'));if ~isempty(Input), obj.Correspondence = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'Kappa'));if ~isempty(Input), obj.Kappa = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'Temp'));if ~isempty(Input), obj.Temp = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'Target'));if ~isempty(Input), obj.Target = varargin{Input+1}; end
          end
        end
    end  
    methods % Special Setting and Getting
        function nr = get.nrPdfP(obj)
            nr = length(obj.PdfP);
        end
        function nr = get.nrSM(obj)
            if isempty(obj.Smeasure), nr = 0; return; end
            nr = obj.Smeasure.nrSM;
        end
        function fields = get.Fields(obj)
            if isempty(obj.Smeasure), fields = {}; return; end
            fields = obj.Smeasure.Fields;
        end
        function obj = set.Fields(obj,field) %#ok<INUSD>
                 % dummy
        end
        function nr = get.nrFields(obj)
            nr = length(obj.Fields);
        end
        function obj = set.Correspondence(obj,cor)
            if isempty(obj.Smeasure), return; end
            obj.Smeasure.Correspondence = cor;
        end       
        function cor = get.Correspondence(obj)
                 if isempty(obj.Smeasure), cor = []; return; end
                 cor = obj.Smeasure.Correspondence;
        end
        function sm = get.Smeasure(obj)
            sm = obj.Smeasure;
            if ~mySuperClass.isH(sm), sm = []; end
        end
        function obj = set.Smeasure(obj,sm)
            if~isempty(obj.Smeasure), delete(obj.Smeasure); end
            obj.Smeasure = sm;
        end
        function nr = get.nrV(obj)
            if isempty(obj.Smeasure),nr= []; return; end
            nr = obj.Smeasure.nrV;
        end
        function target = get.Target(obj)
            if isempty(obj.Smeasure), target = []; return; end
            target = obj.Smeasure.Target;
        end
        function obj = set.Target(obj,target)
            if isempty(obj.Smeasure), return; end
            obj.Smeasure.Target = target;
        end 
        function cp = get.CompleteP(obj)
            cp = obj.CompleteP;
            if~mySuperClass.isH(cp), cp = []; end
        end
        function obj = set.Temp(obj,val)
            if (val<0)||(val>1), error('Temperature must be between 0 and 1'); end
            obj.Temp = val;
        end
    end       
    methods % InterFace functions
        function copy(obj,cobj)
                 cobj.PdfP = obj.PdfP;
                 cobj.CompleteP = obj.CompleteP;% Parents are not cloned!
                 cobj.Evaluation = obj.Evaluation;
                 cobj.Derivative = obj.Derivative;
                 cobj.PdfEvaluation = obj.PdfEvaluation;
                 cobj.Kappa = obj.Kappa;
                 cobj.Temp = obj.Temp;
                 if isempty(obj.Smeasure), return; end
                 cobj.Smeasure = clone(obj.Smeasure);
        end
        function delete(obj)
            if ~isempty(obj.Smeasure),delete(obj.Smeasure); end
        end
        function out = clear(obj)
                 if nargout == 1
                    obj = clone(obj);
                    out = obj;
                 end
                 obj.Evaluation = [];
                 obj.Derivative = [];
                 obj.pdfEvaluation = [];
                 if~isempty(obj.Smeasure), clear(obj.Smeasure); end
        end
        function struc = obj2struc(obj)
            % converting relevant information
            struc = obj2struc@mySuperClass(obj);
            struc.PdfP = obj.PdfP;
            struc.Evaluation = obj.Evaluation;
            struc.Derivative = obj.Derivative;
            struc.pdfEvaluation = obj.pdfEvaluation;
            struc.Kappa = obj.Kappa;
            struc.Temp = obj.Temp;
            struc.Smeasure = [];
            if~isempty(obj.Smeasure), struc.Smeasure = obj2struc(obj.Smeasure); end
        end
        function obj = struc2obj(obj,struc)
                 obj.Type = struc.Type;
                 obj.PdfP = struc.PdfP;
                 obj.Evaluation = struc.Evaluation;
                 obj.Derivative = struc.Derivative;
                 obj.pdfEvaluation = struc.pdfEvaluation;
                 obj.Kappa = struc.Kappa;
                 obj.Temp = struc.Temp;
                 if ~isempty(struc.Smeasure), obj.Smeasure = struc2obj(eval(struc.Smeasure.Type),struc.Smeasure); end
        end
        function [B,sumB] = getB(obj,Tmodel)
            % Getting B coeffs, == LV values in complete process or ones
            if ~isempty(obj.CompleteP)&&(obj.CompleteP.LatentV.nrV==Tmodel.nrV)
               B = obj.CompleteP.LatentV.Value;
               sumB = obj.CompleteP.LatentV.sumV;
               if isempty(B)
                  B = ones(1,Tmodel.nrV);
                  sumB = sum(B);
               end
            else
               B = ones(1,Tmodel.nrV);
               sumB = sum(B);
            end
        end
        function out = updateCorrespondence(obj,Tmodel)
            % updating correspondences for ICP based optimization
            if isempty(obj.Smeasure), error('Cannot update correspondences IP, no Smeasure set'); end
            if nargout == 1
                out = updateCorrespondence(obj.Smeasure,Tmodel);
            else
                updateCorrespondence(obj.Smeasure,Tmodel);
            end           
        end
        function out = initialize(obj,Tmodel)
             if nargout == 1
                obj = clone(obj);
                out = obj;
             end
             if isempty(obj.Smeasure), return; end
             initialize(obj.Smeasure,Tmodel);
             update(obj,Tmodel);
        end
    end
    methods (Abstract = true) %Abstract Interfaces Functions
            out = eval(obj,Tmodel);% evaluation the inlier P
            out = derivate(obj,Tmodel);% derivating the inlier P
            out = update(obj,Tmodel);% updating the parameters of the inlier P
            out = pdf(obj,Tmodel);% evaluation the pdf of the inlier P         
    end
end % classdef

