classdef NextGenBRIMShapeModel < superClass
    properties
        Model = [];% Model to regress against  
        X = [];% Predictors
        XNames = [];% PredictorNames
        ONames = [];% ObservationNames
        RIPX = [];
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
        nrO;% number of observations
        nrY;% number of shape variables
        nrX;% number of predictor variables
    end
    properties
        BRIMModelType = 'PLSR';
        RIPEstMethod = 'UnWeighted';
        ShapeDistance = 'Mahalanobis';
        RIPNormalize = true;
        RIPOutlier = false;
        TrIndex = [];
        InnerFold = 10;
        OuterFold = 12;
        Bootstrap = true;
        MaxIterations = 6;
        StopCorr = 0.98;
        SamplingRuns = 100;
        Show = true;
        CreateMorphs = true;
        Htest = true;
        OuterModels = [];
        BF = [];
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
        function obj = NextGenBRIMShapeModel(varargin)
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
        function out = get.Average(obj)
           if isempty(obj.Model), out = []; return; end
           out = getScan(obj.Model,obj.YAvg);
        end
        function obj = set.X(obj,in)
                 obj.X = in;
                 %initializeRIPX(obj);
        end
    end
    methods % Concrete Basic Interface & Conversion functions
        function obj = initializeRIPX(obj)
            if isempty(obj.X), return; end
            obj.RIPX = nan*zeros(obj.nrO,obj.nrX);
        end
        function out = checkNrObservations(obj)
            if size(obj.X,1)==obj.nrO
                out = true;
            else
                out = false;
            end
        end
        function export2Excell(obj,filename)
            if nargin < 2
               [filename, pathname,filterindex] = uiputfile({ '*.xls','Excell Workbook';...
                                                           '*.xlsx','Excell Workbook';...
                                                             },'Export As',obj.Tag); %#ok<*NASGU>
                if isequal([filename,pathname],[0,0]),return; end
            end
            headers = cell(1,obj.nrX+1);
            RIPheaders = cell(1,obj.nrX+1);
            headers{1} = 'ID';
            RIPheaders{1} = 'ID';
            for i=1:1:obj.nrX
                headers{i+1} = obj.XNames{i};
                RIPheaders{i+1} = ['RIP-' obj.XNames{i}];
            end
            xlswrite(filename,headers,1,'A1');
            xlswrite(filename,RIPheaders,2,'A1');
            xlswrite(filename,obj.ONames',1,'A2');
            xlswrite(filename,obj.ONames',2,'A2');
            xlswrite(filename,obj.X,1,'B2');
            xlswrite(filename,obj.RIPX,2,'B2');
            xlswrite(filename,obj.ONames',3,'A1');
            xlswrite(filename,obj.Model.Tcoeff,3,'B1');  
        end
        function reduceObservations(obj,names)
                 [ind21,ind2] = vlookup(names,obj.ONames);
                 obj.ONames = obj.ONames(ind21);
                 if ~isempty(obj.Model), obj.Model.Tcoeff = obj.Model.Tcoeff(ind21,:);end
                 if ~isempty(obj.X), obj.X = obj.X(ind21,:); end
                 if ~isempty(obj.RIPX), obj.RIPX = obj.RIPX(ind21,:); end
        end
        function out = testSignificance(obj,index,test,t,filename)
           if nargin < 5 
               export = false;
           else
               export = true;
           end
           nrV = length(index);
           switch test
               case 'ROC'
                   out = nan*zeros(nrV,4);
                   for i=1:1:nrV
                       x = obj.X(:,index(i));
                       ripx = obj.RIPX(:,index(i));
                       good = find(~isnan(x));
                       re = x(good);
                       re(find(re==-1)) = 0; %#ok<FNDSB>
                       [tmp,fig] = permROCAnalysis(re,ripx(good),t,obj.XNames{index(i)});
                       out(i,1) = tmp.auc;
                       out(i,2) = tmp.se;
                       out(i,3) = tmp.pauc;
                       out(i,4) = tmp.ppermauc;
                       if export
                          saveas(fig,['ROC_' obj.XNames{index(i)}],'fig');  
                       end
                   end
                   if export
                      headings = cell(1,5);
                      headings{1} = 'Variable';
                      headings{2} = 'Auc';
                      headings{3} = 'se';
                      headings{4} = 'P';
                      headings{5} = 'Perm P';
                      xlswrite(filename,headings,3,'A1'); 
                      xlswrite(filename,obj.XNames(index)',3,'A2');
                      xlswrite(filename,out,3,'B2');
                   end
               case 'COR'
                   out = nan*zeros(nrV,2);
                   for i=1:1:nrV
                       x = obj.X(:,index(i));
                       ripx = obj.RIPX(:,index(i));
                       good = find(~isnan(x));
                       [out(i,1),out(i,2)] = permCorr(x(good),ripx(good),t);
                   end
                   if export
                      headings = cell(1,3);
                      headings{1} = 'Variable';
                      headings{2} = 'Correlation';
                      headings{3} = 'Perm P';
                      xlswrite(filename,headings,2,'A1');
                      xlswrite(filename,obj.XNames(index)',2,'A2');
                      xlswrite(filename,out,2,'B2');
                      close all;
                   end
               case 'ANOVA2'
               case 'ANOVA3'
                   out = nan*zeros(nrV,8);
                   for i=1:1:nrV
                       try
                           x = obj.X(:,index(i));
                           ripx = obj.RIPX(:,index(i));
                           [tmp,fig] = dGPANOVAAnalysis(x,ripx,t,obj.XNames{index(i)});
                           out(i,1) = tmp.ANOVA.FA;out(i,2) = tmp.ANOVA.ppermFA;
                           out(i,3) = tmp.ANOVAPW(3).FA;out(i,4) = tmp.ANOVAPW(3).ppermFA;
                           out(i,5) = tmp.ANOVAPW(1).FA;out(i,6) = tmp.ANOVAPW(1).ppermFA;
                           out(i,7) = tmp.ANOVAPW(2).FA;out(i,8) = tmp.ANOVAPW(2).ppermFA;
                           close all;
                       catch
                       end
                   end
                   if export
                       headers = cell(1,9);
                       headers{1} = 'Variable';headers{2} = 'AA AB BB F';headers{3} = 'AA AB BB P';
                       headers{4} = 'AA BB F';headers{5} = 'AA BB P';headers{6} = 'AA AB F';
                       headers{7} = 'AA AB P';headers{8} = 'BB AB F';headers{9} = 'BB AB P';
                       xlswrite(filename,headers,1,'A1');
                       xlswrite(filename,obj.XNames(index)',1,'A2');
                       xlswrite(filename,out,1,'B2');
                   end       
           end
        end
        function anal = morphAnalysis(obj,index,Cindex,CValues,t,BF)
            Cindex = setdiff(Cindex,index);
            anal = morphPermAnalysis([obj.RIPX(:,Cindex) obj.RIPX(:,index)],...
                                     obj.ShapeLM,obj.Average,BF,CValues,t); 
        end
    end
    methods (Static = true)
       
    end
end

% function [out] = updateRIP(in,M,R2)
%            % in this implementation I take the reference as the origin
%             n2 = size(in,1); % determine input size
%             in = in';
%             n1 = size(M,1);
%             out = nan*zeros(n1,n2);% allocate memory
%             R2n = R2./repmat(sum(R2,2),1,size(R2,2));
%             for j=1:n1
%                 in = in.*repmat(R2n(j,:)',1,n2);
%                 coeff2 = (M(j,:).*R2n(j,:))';
%                 coeff2 = coeff2/norm(coeff2);
%                 out(j,:) = dot(in,repmat(coeff2,1,n2));
%             end
%        end

%                 for i=1:1:n2
%                     %coeff1 = coeffn(:,i)/norm(coeffn(:,i));
%                     % Angle between direction and model direction
%                     %T = coeff1'*coeff2;
%                     %N = sqrt((coeff1'*coeff1)*(coeff2'*coeff2));
%                     %angle = T/N;
%                     % parallell distance (RIP)
%                     %out(j,i) = angle*dist(i);
%                     coeff1 = coeffn(:,i);
%                     out(j,i) = dot(coeff1,coeff2);
%                 end
