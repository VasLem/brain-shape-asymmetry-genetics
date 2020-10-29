classdef BRIMConditionedShapeModel < superClass
    properties
        Model = [];% Model to regress against  
        X = [];% Predictors
        XNames = [];% PredictorNames
        RIPX = [];
        C = [];% Conditioning Variables
        CNames = [];% Conditioning Names
        RIPC = [];
        YType = 'PC';
    end
    properties (Dependent = true)
        Y;% Responses
        YAvg;% average of responses
        YStd;% std of responses
        XAvg;% average of predictors
        XStd;% std of predictors
        RIPXAvg;% average of the RIP 
        RIPXStd;% std of RIP
        CAvg;
        CStd;
        RIPCAvg;
        RIPCStd;
        nrO;% number of observations
        nrY;% number of shape variables
        nrX;% number of predictor variables
        nrC;% number of conditioning variables
    end
    properties
        BRIMModelType = 'PLSR';
        RIPEstimationMethod = 'PLSR Mahalanobis';
        RIPC2Y = [];
        RIPC2RIPX = [];
        RIPX2Y = [];
        TrIndex = [];
        TrainingInnerFold = 10;
        TrainingOuterFold = 2;
        TrainingIterations = 3;
    end
    properties (Dependent = true)
        Average;% average shape
        ShapePC;% shape PCA coefficients
        ShapeLM;% LM configurations shape
        Index;
    end
    properties (Hidden = true, Transient = true)
        HiddenShapeLM = [];
    end
    methods % Constructor
        function obj = BRIMConditionedShapeModel(varargin)
            obj = obj@superClass(varargin{:});         
        end
    end
    methods % Special Setting & Getting
        function out = get.Index(obj)
           if isempty(obj.TrIndex)
              if isempty(obj.Model)
                  out = [];
              else
                  out = (1:obj.Model.n);
              end
           else
               out = obj.TrIndex;
           end  
        end
        function out = get.Model(obj)
           out = obj.Model;
           if ~superClass.isH(out), out = []; end
        end
        function out = get.ShapePC(obj)
            out = obj.Model.Tcoeff(obj.Index,:);
        end
        function out = get.ShapeLM(obj)
            if isempty(obj.HiddenShapeLM)
                tmp = reconstructTraining(obj.Model)';
                out = tmp(obj.Index,:);
                obj.HiddenShapeLM = out;
            else
                out = obj.HiddenShapeLM;
            end
        end
        function obj = set.Model(obj,in)
           obj.Model = in;
           obj.HiddenShapeLM = []; %#ok<*MCSUP>
        end
        function obj = set.TrIndex(obj,in)
           obj.TrIndex = in;
           obj.HiddenShapeLM = [];
        end
        function out = get.Y(obj)
            switch obj.YType
                case 'PC'
                    out = obj.ShapePC;
                case 'LM'
                    out = obj.ShapeLM;
            end
        end
        function out = get.YAvg(obj)
            if isempty(obj.Model), out = []; return; end
            out = mean(obj.Y);
        end
        function out = get.YStd(obj)
            if isempty(obj.Model), out = []; return; end
            out = std(obj.Y);
        end
        function out = get.XAvg(obj)
            if isempty(obj.X), out = []; return; end
            out = nanmean(obj.X);
        end
        function out = get.XStd(obj)
            if isempty(obj.X), out = []; return; end
            out = nanstd(obj.X);
        end
        function out = get.RIPXAvg(obj)
            if isempty(obj.RIPX), out = []; return; end
            out = nanmean(obj.RIPX);
        end
        function out = get.RIPXStd(obj)
            if isempty(obj.RIPX), out = []; return; end
            out = nanstd(obj.RIPX);
        end
        function out = get.CAvg(obj)
            if isempty(obj.C), out = []; return; end
            out = nanmean(obj.C);
        end
        function out = get.CStd(obj)
            if isempty(obj.C), out = []; return; end
            out = nanstd(obj.C);
        end
        function out = get.RIPCAvg(obj)
            if isempty(obj.RIPC), out = []; return; end
            out = nanmean(obj.RIPC);
        end
        function out = get.RIPCStd(obj)
            if isempty(obj.RIPC), out = []; return; end
            out = nanstd(obj.RIPC);
        end
        function out = get.nrO(obj)
                 out = length(obj.Index);
        end
        function out = get.nrY(obj)
            if isempty(obj.Model), out = 0; return; end
            out = size(obj.Y,2);
        end
        function out = get.nrX(obj)
            if isempty(obj.X), out = 0; return; end
            out = size(obj.X,2);
        end
        function out = get.nrC(obj)
            if isempty(obj.C), out = 0; return; end
            out = size(obj.C,2);
        end
        function out = get.Average(obj)
           if isempty(obj.Model), out = []; return; end
           out = getScan(obj.Model,obj.YAvg);
        end
        function obj = set.X(obj,in)
                 obj.X = in;
                 initializeRIPX(obj);
        end
        function obj = set.C(obj,in)
                 obj.C = in;
                 initializeRIPC(obj);
        end
    end
    methods % Concrete Basic Interface & Conversion functions
        function obj = initializeRIPX(obj)
            if isempty(obj.X), return; end
            obj.RIPX = nan*zeros(obj.nrO,obj.nrX);
        end
        function obj = initializeRIPC(obj)
            if isempty(obj.C), return; end
            obj.RIPC = nan*zeros(obj.nrO,obj.nrC);
        end
        function obj = initializeRIP2YM(obj)
           if isempty(obj.RIPX), return; end
           obj.RIP2YM = nan*zeros(obj.nrX,obj.nrY);
           obj.RIP2YMinv = nan*zeros(obj.nrY,obj.nrX);
        end
        function out = checkNrXObservations(obj)
            if size(obj.X,1)==obj.nrO
                out = true;
            else
                out = false;
            end
        end
        function out = checkNrCObservations(obj)
            if size(obj.C,1)==obj.nrO
                out = true;
            else
                out = false;
            end
        end
        function updateRIPC2Y(obj)
                 [IndVar,DepVar] = eliminateNAN(obj.RIPC,obj.Y);
                 [~,~,~,~,M,obj.RIPC2Y.PcVar,~,tmp] = plsregress(IndVar,DepVar,size(IndVar,2));
                 obj.RIPC2Y.M = M;
                 obj.RIPC2Y.Minv = pinv(M);
                 obj.RIPC2Y.YR = tmp.Yresiduals;
                 obj.RIPC2Y.XR = tmp.Xresiduals;

%                  obj.RIPC2Y.M = M(2:end,:);
%                  obj.RIPC2Y.Minv = pinv(obj.RIPC2Y.M);
%                  obj.RIPC2Y.test = pinv(M);
%                  obj.RIPC2Y.YR = tmp.Yresiduals;
%                  obj.RIPC2Y.Intercept = M(1,:);
        end
        function [out,scan] = dRIPC(obj,in,Y,index)
                if nargin<4, index = (1:obj.nrC);end
                M = obj.RIPC2Y.M(2:end,:);
                dY = in*M(index,:);
                out = Y+dY;
                if nargout==1, return; end
                scan = getScan(obj.Model,out);  
        end
       function [out,scan] = changeRIPC(obj,in,RIPC,Y)
                deltaC = in-RIPC;
                out = dRIPC(obj,deltaC,Y);
                if nargout==1, return; end
                scan = getScan(obj.Model,out); 
       end
       function out = RIPCfromY(obj,in)
                if ~(size(in,1)==1), in = in'; end
                out = in*obj.RIPC2Y.Minv;%+obj.RIPCAvg;
                out = out(2:end);
       end
       function [out,scan] = YfromRIPC(obj,in)
                if ~(size(in,1)==1), in = in'; end
                %out = in*obj.RIPC2Y.M+obj.RIPC2Y.Intercept;
                out = [1 in]*obj.RIPC2Y.M;
                if nargout ==1, return; end
                scan = getScan(obj.Model,out);
       end
       function out = RedYfromRIPC(obj,in,RIPC)
                if nargin < 3, RIPC = RIPCfromY(obj,in); end
                %out = in-(RIPC*obj.RIPC2Y.M+obj.RIPC2Y.Intercept);   
                out = in-YfromRIPC(obj,RIPC);
       end
    end
end
%         function updateM(obj,index,Cindex)
%             if nargin < 3, Cindex = [];end
%             if isempty(obj.RIP2YM), initializeRIP2YM(obj); end
%             nr2Condition = length(Cindex);
%             nr2Update = length(index);
%             if ~nr2Condition==0
%                 IndVar = [obj.RIPX(:,Cindex) obj.RIPX(:,index)];
%             else
%                 IndVar = obj.RIPX(:,index);
%             end       
%             UpdateInd = (nr2Condition+1:nr2Condition+nr2Update);
%             DepVar = obj.Y;
%             switch obj.BRIMModelType
%               case 'PLSR'
%                 [IndVar,DepVar] = eliminateNAN(IndVar,DepVar);
%                 [~,~,~,~,M] = plsregress(IndVar,DepVar,size(IndVar,2));
%                 Mpinv = pinv(M(2:end,:));
%                 obj.RIP2YM(index,:) = M(UpdateInd+1,:);
%                 obj.RIP2YMinv(:,index) = Mpinv(:,UpdateInd);
%               otherwise
%             end  
%         end
%         function updateMOneByOne(obj,index,Cindex)
%            updateM(obj,Cindex);
%            for i=1:1:length(index)
%                updateM(obj,index(i),Cindex);
%            end 
%         end
%         function [out,scan] = dRIPX(obj,in,Y,index)
%                 if nargin<4, index = (1:obj.nrX);end
%                 dY = in*obj.RIP2YM(index,:);
%                 out = Y+dY;
%                 if nargout==1, return; end
%                 scan = getScan(obj.Model,out);     
%         end
%        function [out,scan] = changeRIPX(obj,in,RIPX,Y)
%                 deltaX = in-RIPX;
%                 out = dRIPX(obj,deltaX,Y);
%                 if nargout==1, return; end
%                 scan = getScan(obj.Model,out); 
%        end
%        function out = RIPXFromY(obj,in)
%                 if ~(size(in,1)==1), in = in'; end
%                 %out = in*obj.RIP2YMinv+obj.RIPXAvg;
%                 out = zeros(1,obj.nrX);
%                 for x=1:1:obj.nrX
%                     %out(x) = in*obj.RIP2YMinv(:,x)+obj.RIPXAvg(x); 
%                     [~,~,out(x)] = getDistance(obj.Model,[],in,'mahalanobis',obj.RIP2YM(x,:));
%                 end
%        end


%        function updateM(obj)
%             [~,~,~,~,obj.M,obj.PcVar] = plsregress(obj.X,obj.Y,obj.usecmp);
%        end
%        function updateMinv(obj)
%             [~,~,~,~,obj.Minv,obj.PcVar] = plsregress(obj.Y,obj.X,obj.usecmp);
%        end
%        function updateMpinv(obj)
%             if isempty(obj.M), updateM(obj); end
%             obj.Mpinv = pinv(obj.M(2:end,:));
%        end
%        function update(obj)
%           updateM(obj);
%           updateMinv(obj);
%           updateMpinv(obj);
%        end
%        function out = XfromY(obj,in)
%                 if ~(size(in,1)==1), in = in'; end
%                 out = [1 in]*obj.Minv;
%        end
%        function out = XInvfromY(obj,in)
%                 if ~(size(in,1)==1), in = in'; end
%                 out = in*obj.Mpinv+obj.XAvg;
%        end
%        function [out,scan] = YfromX(obj,in)
%                 if ~(size(in,1)==1), in = in'; end
%                 out = [1 in]*obj.M;
%                 if nargout==1, return; end
%                 scan = getScan(obj.Model,out);
%        end
%        function [out,scan] = dX(obj,in,Y,index)
%                 if nargin<4, index = (1:obj.nrX);end
%                 m = obj.M(2:end,:);
%                 dY = in*m(index,:);
%                 out = Y+dY;
%                 if nargout==1, return; end
%                 scan = getScan(obj.Model,out);     
%        end
%        function out = dY(obj,in,X)
%                 dX = in*obj.Minv(2:end,:);
%                 out = X+dX;   
%        end
%        function [out,scan] = changeX(obj,in,X,Y)
%                 deltaX = in-X;
%                 out = dX(obj,deltaX,Y);
%                 if nargout==1, return; end
%                 scan = getScan(obj.Model,out); 
%        end
%        function out = changeY(obj,in,X,Y)
%                 deltaY = in-Y;
%                 out = dY(obj,deltaY,X);
%        end
%        function [outY,outX,scan] = changeSingleX(obj,in,X,Y,index)
%                 dX = X(index)-in;
%                 outY = dX(obj,dX,Y,index);
%                 outX = XfromY(obj,Y);
%                 if nargout==1, return; end
%                 scan = getScan(obj.Model,out);
%        end