classdef shapePLSR < superClass
    properties
        Model = [];% Model to regress against
        TrIndex = [];
        Pred = [];% Predictors
        PredNames = [];% PredictorNames
        nrcmp = [];% nrcmp
        M = [];
        Minv = [];
        PcVar = [];
    end
    properties (Dependent = true)
        Resp;% Responses
        RespAvg;% average of responses
        RespStd;% std of responses
        PredAvg;% average of predictors
        PredStd;% std of predictors
        Index;
        nrO;
        nrResp;
        nrPred;
        usecmp;
    end
    methods % Constructor
        function obj = shapePLSR(varargin)
            obj = obj@superClass(varargin{:});         
        end
    end
    methods % Special Setting & Getting
        function out = get.Index(obj)
           if isempty(obj.TrIndex)
              if obj.nrO==0
                  out = [];
              else
                  out = (1:obj.nrO);
              end
           else
               out = obj.TrIndex;
           end  
        end
        function out = get.usecmp(obj)
           if ~isempty(obj.nrcmp), out = obj.nrcmp; return; end
           out = min(obj.nrResp,obj.nrPred);
        end
        function out = get.Model(obj)
           out = obj.Model;
           if ~superClass.isH(out), out = []; end
        end
        function out = get.Resp(obj)
            out = obj.Model.Tcoeff(obj.Index,:);
        end
        function out = get.RespAvg(obj)
            if isempty(obj.Model), out = []; return; end
            out = mean(obj.Resp);
        end
        function out = get.RespStd(obj)
            if isempty(obj.Model), out = []; return; end
            out = std(obj.Resp);
        end
        function out = get.PredAvg(obj)
            if isempty(obj.Pred), out = []; return; end
            out = mean(obj.Pred);
        end
        function out = get.PredStd(obj)
            if isempty(obj.Pred), out = []; return; end
            out = std(obj.Pred);
        end
        function out = get.nrO(obj)
                 out = lenght(obj.Index);
        end
        function out = get.nrResp(obj)
            if isempty(obj.Model), out = 0; return; end
            out = size(obj.Resp,2);
        end
        function out = get.nrPred(obj)
            if isempty(obj.Pred), out = 0; return; end
            out = size(obj.Pred,2);
        end
    end
    methods % Concrete Basic Interface & Conversion functions
       function out = Coeff2Vec(obj,in)
            % Convert given PCA coeff into a Vector
            if size(in,1)==1, in = in'; end
            if obj.Centering
                out = obj.AvgVec + obj.EigVec*in;
            else
                out = obj.EigVec*in;
            end
       end
    end
end