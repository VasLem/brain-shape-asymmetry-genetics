classdef DNA2FaceModel < superClass
    properties
        Model = [];% Model to regress against
        TrIndex = [];
        Sex = [];%Conditioning Variables
        RIPS = [];
        Ancestry = [];
        RIPA = [];
        Genes = [];
        RIPG = [];
        GenesNames =  []% Gen Names 
        nrcmp = [];% nrcmp
        InnerFold = 20;
        OuterFold = 2;
%         M = [];
%         Minv = [];
%         Mpinv = [];
%         PcVar = [];
    end
    properties (Dependent = true)
        Average;
        ShapePC;
        ShapeLM;
        Y;% Responses
        YAvg;% average of responses
        YStd;% std of responses
        X;% Predictors
        XNames;% PredictorNames
        XAvg;% average of predictors
        XStd;% std of predictors
        Index;
        nrO;
        nrY;
        nrX;
        usecmp;
    end
    methods % Constructor
        function obj = DNA2FaceModel(varargin)
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
        function out = get.X(obj)
           out = [obj.Sex obj.Ancestry obj.Genes];
        end
        function out = get.XNames(obj)
           out = ['Sex' 'Ancestry' obj.GenesNames]; 
        end
        function out = get.usecmp(obj)
           if ~isempty(obj.nrcmp), out = obj.nrcmp; return; end
           out = min(obj.nrY,obj.nrX);
        end
        function out = get.Model(obj)
           out = obj.Model;
           if ~superClass.isH(out), out = []; end
        end
        function out = get.Y(obj)
            out = obj.Model.Tcoeff(obj.Index,:);
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
            out = mean(obj.X);
        end
        function out = get.XStd(obj)
            if isempty(obj.X), out = []; return; end
            out = std(obj.X);
        end
        function out = get.nrO(obj)
                 out = lenght(obj.Index);
        end
        function out = get.nrY(obj)
            if isempty(obj.Model), out = 0; return; end
            out = size(obj.Y,2);
        end
        function out = get.nrX(obj)
            if isempty(obj.X), out = 0; return; end
            out = size(obj.X,2);
        end
        function out = get.Average(obj)
           if isempty(obj.Model), out = []; return; end
           out = getScan(obj.Model,obj.YAvg);
        end
        function out = get.ShapePC(obj)
           out =  obj.Model.Tcoeff(obj.Index,:);
        end
    end
    methods % Concrete Basic Interface & Conversion functions
       function updateM(obj)
            [~,~,~,~,obj.M,obj.PcVar] = plsregress(obj.X,obj.Y,obj.usecmp);
       end
       function updateMinv(obj)
            [~,~,~,~,obj.Minv,obj.PcVar] = plsregress(obj.Y,obj.X,obj.usecmp);
       end
       function updateMpinv(obj)
            if isempty(obj.M), updateM(obj); end
            obj.Mpinv = pinv(obj.M(2:end,:));
       end
       function update(obj)
          updateM(obj);
          updateMinv(obj);
          updateMpinv(obj);
       end
       function out = XfromY(obj,in)
                if ~(size(in,1)==1), in = in'; end
                out = [1 in]*obj.Minv;
       end
       function out = XInvfromY(obj,in)
                if ~(size(in,1)==1), in = in'; end
                out = in*obj.Mpinv+obj.XAvg;
       end
       function [out,scan] = YfromX(obj,in)
                if ~(size(in,1)==1), in = in'; end
                out = [1 in]*obj.M;
                if nargout==1, return; end
                scan = getScan(obj.Model,out);
       end
       function [out,scan] = dX(obj,in,Y,index)
                if nargin<4, index = (1:obj.nrX);end
                m = obj.M(2:end,:);
                dY = in*m(index,:);
                out = Y+dY;
                if nargout==1, return; end
                scan = getScan(obj.Model,out);     
       end
       function out = dY(obj,in,X)
                dX = in*obj.Minv(2:end,:);
                out = X+dX;   
       end
       function [out,scan] = changeX(obj,in,X,Y)
                deltaX = in-X;
                out = dX(obj,deltaX,Y);
                if nargout==1, return; end
                scan = getScan(obj.Model,out); 
       end
       function out = changeY(obj,in,X,Y)
                deltaY = in-Y;
                out = dY(obj,deltaY,X);
       end
       function [outY,outX,scan] = changeSingleX(obj,in,X,Y,index)
                dX = X(index)-in;
                outY = dX(obj,dX,Y,index);
                outX = XfromY(obj,Y);
                if nargout==1, return; end
                scan = getScan(obj.Model,out);
       end
    end
end