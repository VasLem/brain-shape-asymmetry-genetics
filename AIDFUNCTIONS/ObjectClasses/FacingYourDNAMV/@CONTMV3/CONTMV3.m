classdef CONTMV3 < CATMV3
    properties
    end
    properties (Dependent = true)
    end
    properties (Hidden = true) 
    end
    properties(Hidden = true, Dependent = true)
        XREG;
        PosClass;
    end
    methods % CONSTRUCTOR
        function obj = CONTMV3(varargin)
            obj = obj@CATMV3(varargin{:});         
        end
    end
    methods % GENERAL SETTING/GETTING
        function out = get.XREG(obj)
            out = obj.X;
        end
        function out = get.PosClass(obj) %#ok<*MANU>
           out = 1; 
        end
    end
    methods % PREDICTOR
    end
    methods % BIOMETRICS
    end
    methods % EVOMORPH
    end
    methods % INTERFACING
    end
    methods (Static = true)
    end
end


% function [NX,COST,W,ind,Tind,Find] = prepPredictorTraining(obj,X,TargetX)
%             obj.TargX = TargetX;
%             NX = nan*zeros(size(X));
%             W = nan*zeros(size(X));
%             ind = find(~isnan(X));
%             
%             range = max(X)-min(X);
%             range = range/obj.Bins;
%             obj.RangeX = range;
%             %medX = median(X(ind));
%             %madX = mad(X(ind));
%             %range = mad(X(ind));
%            % select closest examples
%             dif = abs(TargetX-X(ind));
%             Tind2 = find(dif<=range);
%             Tind = ind(Tind2);
%             %[~,ind2] = sort(dif,'ascend');
%             %Tind2 = ind2(1:obj.KN);
%             %Tind = ind(Tind2);
%             NX(Tind) = obj.PosClass;
%             tmp = dif(Tind2);
%             if max(tmp)==0
%                W(Tind) = 1;
%             else
%                tmp = tmp-min(tmp);tmp = tmp/max(tmp);tmp = 1-tmp;
%                %mval = min(tmp(tmp>0));tmp(tmp==0) = mval;
%                W(Tind) = tmp;
%             end
%            % opposite examples
%             Find2 = find(dif>range);
%             %Find2 = ind2(end-obj.KN+1:end);
%             Find = ind(Find2);
%             NX(Find) = -1*obj.PosClass;
%             tmp = dif(Find2);
%             tmp = tmp-min(tmp);tmp = tmp/max(tmp);
%             %mval = min(tmp(tmp>0));tmp(tmp==0) = mval;
%             W(Find) = tmp;
%             ind = union(Tind,Find);
%             W(ind) = W(ind)+1;
%            % Defining COST 
%             COST.ClassNames = [obj.PosClass -1*obj.PosClass];
%             COST.ClassificationCosts = zeros(2,2);
%             %COST.ClassificationCosts(1,2) = 1;
%             %COST.ClassificationCosts(2,1) = 1;
%             COST.ClassificationCosts(1,2) = 1-(length(Tind)/length(ind));
%             COST.ClassificationCosts(2,1) = 1-(length(Find)/length(ind));
%             W = nan*zeros(size(X));
%             W(Tind) = 1-(length(Tind)/length(ind));
%             W(Find) = 1-(length(Find)/length(ind));
%         end

% function [MS,MH] = match2Seq(obj,TargX,FS)
%             if isempty(obj.SeqX), error('Sequence not trained'); end
%             n = length(obj.SeqX);
%             diff = zeros(1,length(obj.SeqX));
%             in = zeros(1,length(obj.SeqX));
%             nT = size(FS,1);
%             MS = zeros(nT,n);
%             MH = zeros(nT,n);
%             for i=1:n
%                 diff(i) = abs(obj.SeqX{i}.TargX-TargX);
%                 if diff(i)<=obj.SeqX{i}.RangeX
%                    in(i) = 1;
%                 else
%                    in(i) = -1;%continue;
%                 end
%                 [PX,PXC] = predict(obj.SeqX{i},FS);
%                 [MS(:,i),MH(:,i)] = Pred2Match(obj.SeqX{i},in(i),PX,PXC);
%             end
%             %good = find(in==1);
%             good = 1:n;
%             W = 1./diff(good);
%             sumW = sum(W);
%             MS = MS(:,good);
%             MS = sum(repmat(W,nT,1).*MS,2)/sumW;
%             MH = MH(:,good);
%             MH = sum(repmat(W,nT,1).*MH,2)/sumW;  
%         end
