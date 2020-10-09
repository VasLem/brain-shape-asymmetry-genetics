classdef PLSRSuper < superClass
    properties (Abstract = true)
        Y;% Model to regress against
        X;% Predictors
    end
    properties (Dependent = true)
        YAvg;% average of responses
        YStd;% std of responses
        XAvg;% average of predictors
        XStd;% std of predictors
        nrO;
        nrY;
        nrX;
        usecmp;
    end
    properties
        M = [];
        Minv = [];
        PcVar = []; 
        nrcmp = [];% nrcmp
        YResiduals = [];
        XResiduals = [];
    end
    methods % Constructor
        function obj = PLSRSuper(varargin)
            obj = obj@superClass(varargin{:});         
        end
    end
    methods % Special Setting & Getting
        function out = get.usecmp(obj)
           if ~isempty(obj.nrcmp), out = obj.nrcmp; return; end
           out = min(obj.nrY,obj.nrX);
        end
        function out = get.YAvg(obj)
            if isempty(obj.Y), out = []; return; end
            out = nanmean(obj.Y);
        end
        function out = get.YStd(obj)
            if isempty(obj.Y), out = []; return; end
            out = nanstd(obj.Y);
        end
        function out = get.XAvg(obj)
            if isempty(obj.X), out = []; return; end
            out = nanmean(obj.X);
        end
        function out = get.XStd(obj)
            if isempty(obj.X), out = []; return; end
            out = nanstd(obj.X);
        end
        function out = get.nrO(obj)
                 out = size(obj.Y,1);
        end
        function out = get.nrY(obj)
            out = size(obj.Y,2);
        end
        function out = get.nrX(obj)
            if isempty(obj.X), out = 0; return; end
            out = size(obj.X,2);
        end
    end
    methods % Concrete Basic Interface & Conversion functions
       function out = checkNrObservations(obj)
            if size(obj.X,1)==obj.nrO
                out = true;
            else
                out = false;
            end
        end 
       function update(obj)
            if ~checkNrObservations(obj), error('different amount of observations between X and Y'); end
            [IndVar,DepVar,ind1] = eliminateNAN(obj.X,obj.Y);
            [DepVar,IndVar,ind2] = eliminateNAN(DepVar,IndVar);
            [~,~,~,~,obj.M,obj.PcVar,~,tmp] = plsregress(IndVar,DepVar,obj.usecmp);
            obj.Minv = pinv(obj.M);
            obj.XResiduals = nan*zeros(size(obj.X));
            obj.YResiduals = nan*zeros(size(obj.Y));
            obj.XResiduals(ind1(ind2),:) = tmp.Xresiduals;
            obj.YResiduals(ind1(ind2),:) = tmp.Yresiduals;
       end
       function out = dX(obj,in,Y,index)
                if nargin<4, index = (1:obj.nrX);end
                m = obj.M(2:end,:);
                dY = in*m(index,:);
                out = Y+dY;    
       end
       function out = changeX(obj,in,X,Y)
                deltaX = in-X;
                out = dX(obj,deltaX,Y,1:obj.nrX);
       end
       function out = XfromY(obj,in)
                if ~(size(in,1)==1), in = in'; end
                out = in*obj.Minv;
                out = out(2:end);
       end
       function out = YfromX(obj,in)
                if ~(size(in,1)==1), in = in'; end
                out = [1 in]*obj.M;
       end
       function out = YResfromX(obj,in,X)
           if nargin < 3, X = XfromY(obj,in); end
           out = in-YfromX(obj,X);
       end
       function out = YRes(obj,in)
           x = XfromY(obj,in);
           out = in-YfromX(obj,x);
       end
    end
end