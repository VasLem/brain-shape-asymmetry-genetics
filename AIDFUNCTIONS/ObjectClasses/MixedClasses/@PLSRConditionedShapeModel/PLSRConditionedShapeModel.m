classdef PLSRConditionedShapeModel < PLSRSuper
    properties (Dependent = true)
       Model; % stored in CYModel
    end
    properties
        CYModel = PLSRShapeModel;
        CXModel = [];
        XYModel = [];
        X = [];
        XNames = {};
    end
    properties (Dependent = true)
        %X;% stored in CXModel
        %XNames;% stored in XYModel
        Y;% stored in CYModel
        C;% stored in CYModel
        CNames;% Conditioning variable names, stored in CYModel
        CAvg;% average of conditioning
        CStd;% Std of conditioning
        nrC;% nr of conditioning variables
        YCResiduals;% Y residuals after conditioning on the CondModel
        XCResiduals;% X residuals after regressing C onto X.
        Average;% stored in CYModel
    end
    methods % Constructor
        function obj = PLSRConditionedShapeModel(varargin)
            obj = obj@PLSRSuper(varargin{:});
            obj.CXModel = PLSRC2X;
            obj.XYModel = PLSRReduced;
            obj.CXModel.Parent = obj;
            obj.XYModel.Parent = obj;
        end
    end
    methods % Special Setting & Getting
        function out = get.CYModel(obj)
           out = obj.CYModel;
           if ~superClass.isH(out)
               new = PLSRShapeModel;
               obj.CYModel = new;
               out = new;
           end
        end
        function out = get.CXModel(obj)
           out = obj.CXModel;
           if ~superClass.isH(out)
               new = PLSRC2X;
               new.Parent = obj;
               obj.CXModel = new;
               out = new;
           end
        end
        function out = get.XYModel(obj)
            out = obj.XYModel;
            if ~superClass.isH(out)
               new = PLSRReduced;
               new.Parent = obj;
               obj.XYModel = new;
               out = new;
           end
        end
        function out = get.Model(obj)
           if isempty(obj.CYModel), out = []; return; end
           out = obj.CYModel.Model;
           if ~superClass.isH(out), out = []; end
        end
        function obj = set.Model(obj,in)
                 obj.CYModel.Model = in;
        end
        function out = get.Y(obj)
            out = obj.CYModel.Y;
        end
        function out = get.Average(obj)
           out = obj.CYModel.Average;
        end
        function obj = set.C(obj,in)
            obj.CYModel.X = in;
        end
        function out = get.C(obj)
            out = obj.CYModel.X;
        end
        function obj = set.CNames(obj,in)
            obj.CYModel.XNames = in;
        end
        function out = get.CNames(obj)
            out = obj.CYModel.XNames;
        end
        function out = get.CAvg(obj)
           out = obj.CYModel.XAvg; 
        end
        function out = get.CStd(obj)
           out = obj.CYModel.XStd;
        end
        function out = get.nrC(obj)
            out = obj.CYModel.nrX;
        end
        function out = get.YCResiduals(obj)
            out = obj.CYModel.YResiduals;
        end
        function out = get.XCResiduals(obj)
            out = obj.CXModel.YResiduals;
        end
    end
    methods % Concrete Basic Interface & Conversion functions
        function out = checkNrCObservations(obj)
            out = checkNrObservations(obj.CYModel);
        end
        function update(obj)
           updateCYModel(obj);
           updateCXModel(obj);
           updateXYModel(obj);
        end
        function updateX(obj)
            updateCXModel(obj);
            updateXYModel(obj);
        end
        function updateCYModel(obj)           
            update(obj.CYModel);
        end
        function updateCXModel(obj)
           update(obj.CXModel); 
        end
        function updateXYModel(obj)
            update(obj.XYModel);
        end
        function [out,scan] = dC(obj,in,Y,index)
                if nargin<4, index = (1:obj.nrC);end
                [out,scan] = dX(obj.CYModel,in,Y,index); 
        end
        function [out,scan] = changeCY(obj,in,C,Y)
                [out,scan] = changeX(obj.CYModel,in,C,Y);
        end
        function [out] = changeCX(obj,in,C,X)
            out = changeX(obj.CXModel,in,C,X);
        end
        function [outY,outX,scan] = changeC(obj,in,C,Y,X)
                [outY,scan] = changeX(obj.CYModel,in,C,Y);
                [outX] = changeX(obj.CXModel,in,C,X);
        end
        function out = CfromY(obj,in)
                out = XfromY(obj.CYModel,in);
        end
        function [out,scan] = YfromC(obj,in)
                [out,scan] = YfromX(obj.CYModel,in);
        end
        function out = XfromC(obj,in)
            out = YfromX(obj.CXModel,in);
        end
        function out = CfromX(obj,in)
           out = XfromY(obj.CXModel,in); 
        end
        function out = XfromYC(obj,Y,C)
           yres =  YResfromC(obj,Y,C);
           xres = XfromY(obj.XYModel,yres);
           xbase = XfromC(obj,C);
           out = xbase + xres;
        end
        function [out,scan] = YfromXC(obj,X,C)
            ybase = YfromC(obj,C);
            xres = XResfromC(obj,X,C);
            yres = YfromX(obj.XYModel,xres);
            out = ybase + yres;
            if nargout==1, return; end
            scan = getScan(obj.Model,out);
        end
        function out = YResfromC(obj,in,C)
           out = YResfromX(obj.CYModel,in,C);
        end
        function out = XResfromC(obj,in,C)
            out = YResfromX(obj.CXModel,in,C);
        end
        function [out,scan] = changeX(obj,in,X,Y,C)
            % get Y residual from C
                ybase = YfromC(obj,C);
                ycurres = Y-ybase;
            % get X residual current and in from C
                xres = XResfromC(obj,in,C);
                xcurres = XResfromC(obj,X,C);
            % change X current to X in in reduced model space and create out  
                yres = changeX(obj.XYModel,xres,xcurres,ycurres);
                out = ybase + yres;
            if nargout==1, return; end
            scan = getScan(obj.Model,out);
        end
        function reconnectChildren(obj)
            obj.CXModel.Parent = obj;
            obj.XYModel.Parent = obj;
        end
    end
end




%         function obj = set.XNames(obj,in)
%            obj.XYModel.XNames = in; 
%         end
%         function out = get.XNames(obj)
%             out = obj.XYModel.XNames;
%         end

%         function obj = set.X(obj,in)
%             obj.CXModel.Y = in;
%         end
%         function out = get.X(obj)
%             out = obj.CXModel.Y;
%         end


    %        function out = checkNrXObservations(obj)
%             if size(obj.X,1)==obj.nrO
%                 out = true;
%             else
%                 out = false;
%             end
%         end 
%        function update(obj)
%             if ~checkNrObservations(obj), error('different amount of observations between X and Y'); end
%             [IndVar,DepVar] = eliminateNAN(obj.X,obj.Y);
%             [~,~,~,~,obj.M,obj.PcVar,~,tmp] = plsregress(IndVar,DepVar,obj.usecmp);
%             obj.Minv = pinv(obj.M);
%             obj.XResiduals = tmp.Xresiduals;
%             obj.YResiduals = tmp.Yresiduals;
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



%        function [outY,outX,scan] = changeSingleX(obj,in,X,Y,index)
%                 dX = X(index)-in;
%                 outY = dX(obj,dX,Y,index);
%                 outX = XfromY(obj,Y);
%                 if nargout==1, return; end
%                 scan = getScan(obj.Model,out);
%        end
