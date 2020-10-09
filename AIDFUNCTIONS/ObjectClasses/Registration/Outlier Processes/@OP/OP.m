classdef OP < mySuperClass
    % This is the abstract interface class for Transformation models
    properties 
        CompleteP = []; % link to Parent, being the complete Process
        Evaluation = [];% Outlier Process negative log evaluations over all floating surface points
        Derivative = [];% Outlier Process gradients in floating surface points
        PdfEvaluation = [];% Outlier Process pdf evaluation of floating surface points
    end
    properties (Abstract = true)
        PdfP;% Outlier Process parameters
    end
    properties (Dependent = true)
        nrPdfP;% number of parameters
        nrV;% number of point evaluations
    end
    methods %Constructor
        function obj = OP(varargin)
            obj = obj@mySuperClass(varargin{:});
          if nargin > 0
             Input = find(strcmp(varargin, 'PdfP'));if ~isempty(Input), obj.PdfP = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'CompleteP'));if ~isempty(Input), obj.CompleteP = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'Evaluation'));if ~isempty(Input), obj.Evaluation = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'Derivative'));if ~isempty(Input), obj.Derivative = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'PdfEvaluation'));if ~isempty(Input), obj.PdfEvaluation = varargin{Input+1}; end
          end
        end
    end  
    methods % Special Setting and Getting
        function nr = get.nrPdfP(obj)
            nr = length(obj.PdfP);
        end
        function cp = get.CompleteP(obj)
            cp = obj.CompleteP;
            if~mySuperClass.isH(cp), cp = []; end
        end
        function nr = get.nrV(obj)
            nr = length(obj.PdfEvaluation);
        end
    end       
    methods % InterFace functions 
        function copy(obj,cobj)
                 copy@mySuperClass(obj,cobj)
                 cobj.PdfP = obj.PdfP;
                 cobj.CompleteP = obj.CompleteP;% Parents are not cloned!
                 cobj.Evaluation = obj.Evaluation;
                 cobj.Derivative = obj.Derivative;
                 cobj.PdfEvaluation = obj.PdfEvaluation;         
        end
        function varargout = clear(obj)
                 if nargout == 1
                    obj = clone(obj);
                    varargout{1} = obj;
                 end
                 obj.Evaluation = [];
                 obj.Derivative = [];
                 obj.PdfEvaluation = [];
        end
        function struc = obj2struc(obj)
            % converting relevant information
            struc = obj2struc@mySuperClass(obj);
            struc.PdfP = obj.PdfP;
            struc.Evaluation = obj.Evaluation;
            struc.Derivative = obj.Derivative;
            struc.PdfEvaluation = obj.PdfEvaluation;
        end
        function obj = struc2obj(obj,struc)
                 obj.PdfP = struc.PdfP;
                 obj.Evaluation = struc.Evaluation;
                 obj.Derivative = struc.Derivative;
                 obj.PdfEvaluation = struc.PdfEvaluation;
        end
        function [B,sumB] = getB(obj,Tmodel)
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
        function out = initialize(obj,Tmodel)
             if nargout == 1
                obj = clone(obj);
                out = obj;
             end
             update(obj,Tmodel);
        end
    end
    methods (Abstract = true) %Abstract Interfaces Functions
            out = eval(obj,Tmodel);
            out = derivate(obj,Tmodel);
            out = update(obj,Tmodel);
            out = pdf(obj,Tmodel);
    end
end % classdef

