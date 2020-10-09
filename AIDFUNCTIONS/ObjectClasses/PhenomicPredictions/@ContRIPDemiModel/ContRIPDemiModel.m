classdef ContRIPDemiModel < RIPDemiModel
    properties
        Method = 'polyfit';
        TrP;
        TrS;
        D = 1;
    end
    methods % Constructor
        function obj = ContRIPDemiModel(varargin)
            obj = obj@RIPDemiModel(varargin{:});         
            obj.XType = 'Continous';
            obj.Level = 'Genomic';% most of the time this is true so set this
        end
    end
    methods % Special Setting & Getting
    end
    methods % General Interface Functions
        function out = getScore(obj,in,val)
           out = getRIP(obj,in);
           [out,E] = polyval(obj.TrP,out',obj.TrS);
           out = abs(out-val)./E;
        end
        function trainModel(obj)
            createMEstRIP(obj)
            % optimal threshold is obtained using ROC analysis
            [X,RIP] = eliminateNAN(obj.TrX,obj.TrRIP);
            [obj.TrP,obj.TrS] = polyfit(RIP,X,obj.D);
            % Test the model accuracy
            getAcc(obj);
        end
        function deployModel(obj)
           trainModel(obj);
           obj.TrY = [];
           %obj.TrRIP = [];
           %obj.TrX = [];
           %obj.EstRIP = [];
        end
        function getOutSampleAcc(obj,Y,X)
            index = ~isnan(X);
            X = X(index);
            Y = Y(index,:);
            rip = getRIP(obj,Y');
            [out,~] = polyval(obj.TrP,rip',obj.TrS);
            tmp = corrcoef(X,out);
            obj.OutSampleAcc = tmp(1,2);
        end
    end
end