classdef PLSRShapeModel < PLSRSuper
    properties
        Model = [];% Model to regress against
        X = [];% Predictors
        XNames = [];% PredictorNames
    end
    properties (Dependent = true)
        Y;% Responses
        Average;
        %ShapePC;% shape PCA coefficients
        %ShapeLM;% LM configurations shape
    end
    methods % Constructor
        function obj = PLSRShapeModel(varargin)
            obj = obj@PLSRSuper(varargin{:});         
        end
    end
    methods % Special Setting & Getting
        function out = get.Model(obj)
           out = obj.Model;
           if ~superClass.isH(out), out = []; end
        end
        function out = get.Y(obj)
            if isempty(obj.Model), out = []; return; end
            out = obj.Model.Tcoeff;
        end
        function out = get.Average(obj)
           if isempty(obj.Model), out = []; return; end
           out = getScan(obj.Model,obj.YAvg);
        end
    end
    methods % Concrete Basic Interface & Conversion functions
       function [out,scan] = dX(obj,in,Y,index)
                out = dX@PLSRSuper(obj,in,Y,index);
                if nargout==1, return; end
                scan = getScan(obj.Model,out);     
       end
       function [out,scan] = changeX(obj,in,X,Y)
                out = changeX@PLSRSuper(obj,in,X,Y);
                if nargout==1, return; end
                scan = getScan(obj.Model,out); 
       end
       function [out,scan] = YfromX(obj,in)
                out = YfromX@PLSRSuper(obj,in);
                if nargout==1, return; end
                scan = getScan(obj.Model,out);
       end      
    end
end


%        function [outY,outX,scan] = changeSingleX(obj,in,X,Y,index)
%                 dX = X(index)-in;
%                 outY = dX(obj,dX,Y,index);
%                 outX = XfromY(obj,Y);
%                 if nargout==1, return; end
%                 scan = getScan(obj.Model,out);
%        end


% classdef PLSRShapeModel < superClass
%     properties
%         Model = [];% Model to regress against
%         X = [];% Predictors
%         XNames = [];% PredictorNames
%     end
%     properties (Dependent = true)
%         Y;% Responses
%         YAvg;% average of responses
%         YStd;% std of responses
%         XAvg;% average of predictors
%         XStd;% std of predictors
%         Index;
%         nrO;
%         nrY;
%         nrX;
%         usecmp;
%         Average;
%         %ShapePC;% shape PCA coefficients
%         %ShapeLM;% LM configurations shape
%     end
%     properties
%         M = [];
%         Minv = [];
%         PcVar = []; 
%         nrcmp = [];% nrcmp
%         TrIndex = [];
%         YResiduals = [];
%         XResiduals = [];
%     end
%     methods % Constructor
%         function obj = PLSRShapeModel(varargin)
%             obj = obj@superClass(varargin{:});         
%         end
%     end
%     methods % Special Setting & Getting
%         function out = get.Index(obj)
%            if isempty(obj.TrIndex)
%               if isempty(obj.Model)
%                   out = [];
%               else
%                   out = (1:obj.Model.n);
%               end
%            else
%                out = obj.TrIndex;
%            end  
%         end
%         function out = get.usecmp(obj)
%            if ~isempty(obj.nrcmp), out = obj.nrcmp; return; end
%            out = min(obj.nrY,obj.nrX);
%         end
%         function out = get.Model(obj)
%            out = obj.Model;
%            if ~superClass.isH(out), out = []; end
%         end
%         function out = get.Y(obj)
%             out = obj.Model.Tcoeff(obj.Index,:);
%         end
%         function out = get.YAvg(obj)
%             if isempty(obj.Model), out = []; return; end
%             out = mean(obj.Y);
%         end
%         function out = get.YStd(obj)
%             if isempty(obj.Model), out = []; return; end
%             out = std(obj.Y);
%         end
%         function out = get.XAvg(obj)
%             if isempty(obj.X), out = []; return; end
%             out = nanmean(obj.X);
%         end
%         function out = get.XStd(obj)
%             if isempty(obj.X), out = []; return; end
%             out = nanstd(obj.X);
%         end
%         function out = get.nrO(obj)
%                  out = length(obj.Index);
%         end
%         function out = get.nrY(obj)
%             if isempty(obj.Model), out = 0; return; end
%             out = size(obj.Y,2);
%         end
%         function out = get.nrX(obj)
%             if isempty(obj.X), out = 0; return; end
%             out = size(obj.X,2);
%         end
%         function out = get.Average(obj)
%            if isempty(obj.Model), out = []; return; end
%            out = getScan(obj.Model,obj.YAvg);
%         end
%     end
%     methods % Concrete Basic Interface & Conversion functions
%        function out = checkNrObservations(obj)
%             if size(obj.X,1)==obj.nrO
%                 out = true;
%             else
%                 out = false;
%             end
%         end 
%        function update(obj)
%             if ~checkNrObservations(obj), error('different amount of observations between X and Y'); end
% %             [IndVar,DepVar] = eliminateNAN(obj.X,obj.Y);
% %             [~,~,~,~,obj.M,obj.PcVar,~,tmp] = plsregress(IndVar,DepVar,obj.usecmp);
% %             obj.Minv = pinv(obj.M);
% %             obj.XResiduals = tmp.Xresiduals;
% %             obj.YResiduals = tmp.Yresiduals;
%             [IndVar,DepVar,ind1] = eliminateNAN(obj.X,obj.Y);
%             [DepVar,IndVar,ind2] = eliminateNAN(DepVar,IndVar);
%             [~,~,~,~,obj.M,obj.PcVar,~,tmp] = plsregress(IndVar,DepVar,obj.usecmp);
%             obj.Minv = pinv(obj.M);
%             obj.XResiduals = nan*zeros(size(obj.X));
%             obj.YResiduals = nan*zeros(size(obj.Y));
%             obj.XResiduals(ind1(ind2),:) = tmp.Xresiduals;
%             obj.YResiduals(ind1(ind2),:) = tmp.Yresiduals;
%        end
%        function [out,scan] = dX(obj,in,Y,index)
%                 if nargin<4, index = (1:obj.nrX);end
%                 m = obj.M(2:end,:);
%                 dY = in*m(index,:);
%                 out = Y+dY;
%                 if nargout==1, return; end
%                 scan = getScan(obj.Model,out);     
%        end
%        function [out,scan] = changeX(obj,in,X,Y)
%                 deltaX = in-X;
%                 out = dX(obj,deltaX,Y);
%                 if nargout==1, return; end
%                 scan = getScan(obj.Model,out); 
%        end
%        function out = XfromY(obj,in)
%                 if ~(size(in,1)==1), in = in'; end
%                 out = in*obj.Minv;
%                 out = out(2:end);
%        end
%        function [out,scan] = YfromX(obj,in)
%                 if ~(size(in,1)==1), in = in'; end
%                 out = [1 in]*obj.M;
%                 if nargout==1, return; end
%                 scan = getScan(obj.Model,out);
%        end
%        function out = YResfromX(obj,in,X)
%            if nargin < 3, X = XfromY(obj,in); end
%            out = in-YfromX(obj,X);
%        end       
%     end
% end
% 
% 
% %        function [outY,outX,scan] = changeSingleX(obj,in,X,Y,index)
% %                 dX = X(index)-in;
% %                 outY = dX(obj,dX,Y,index);
% %                 outX = XfromY(obj,Y);
% %                 if nargout==1, return; end
% %                 scan = getScan(obj.Model,out);
% %        end
