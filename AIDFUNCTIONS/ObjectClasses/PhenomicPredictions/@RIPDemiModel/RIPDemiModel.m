classdef RIPDemiModel < demiModel
    properties
        TrRIP;% Given RIP values of training faces
        TrC;% Covariates used during regression
        EstRIP;% Estimated RIP values of training faces
        YEigStd;% Eigenvalues of the shape/appearance model used to train the model
        YAvg;% Average coeff of the shape/appearance model, to use as a reference int estimating RIP values
        TrM;% The resulting regression matrix from the training data.
        OutSampleAcc;% Testing the DemiModel on a sample outside the training
        Flipped;
    end
    properties (Dependent = true)
       nrC;% number of covariates
       TrLOOCor;
       XTrCor;
       XLOOCor;
    end
    methods % Constructor
        function obj = RIPDemiModel(varargin)
            obj = obj@demiModel(varargin{:}); 
        end
    end
    methods % Special Setting & Getting
        function out = get.nrC(obj)
            if isempty(obj.TrC), out = 0; return; end
            out = size(obj.TrC,2);
        end
        function out = get.TrLOOCor(obj)
           if isempty(obj.TrRIP); out = []; return; end
           if isempty(obj.EstRIP); out = []; return; end
           out = permCorr(obj.TrRIP,obj.EstRIP,0);
        end
        function out = get.XTrCor(obj)
           if isempty(obj.TrRIP); out = []; return; end
           if isempty(obj.TrX); out = []; return; end
           index = ~isnan(obj.TrX);
           out = permCorr(obj.TrRIP(index),obj.TrX(index),0);
        end
        function out = get.XLOOCor(obj)
           if isempty(obj.EstRIP); out = []; return; end
           if isempty(obj.TrX); out = []; return; end
           index = ~isnan(obj.TrX);
           out = permCorr(obj.EstRIP(index),obj.TrX(index),0);
        end
    end
    methods % General Interface Functions
        function out = getRIP(obj,in)% fast implementation, can handle a whole population at once
            if isempty(obj.TrM), trainModel(obj); end
            out = RIPDemiModel.getRIPStatic(obj.YAvg,obj.YEigStd,in,obj.TrM);
        end
        function createMEstRIP(obj)          
            if obj.nrC==0
               RIP = obj.TrRIP;
            else
               RIP = [obj.TrC obj.TrRIP];
               disp('With TrC');
            end
            [RIP,Y] = eliminateNAN(RIP,obj.TrY);
            [~,~,~,~,obj.TrM] = plsregress(RIP,Y,size(RIP,2));       
            TrY = obj.TrY;
            EstRIP = zeros(1,obj.nrTr); %#ok<*PROP>
            Ref = obj.YAvg;
            std = obj.YEigStd;
            nrTr = obj.nrTr;
            %tic;
            parfor i=1:nrTr
               build = setdiff((1:nrTr),i); 
               [forRIP,Y] = eliminateNAN(RIP(build,:),TrY(build,:));
               [~,~,~,~,M] = plsregress(forRIP,Y,size(forRIP,2));
               out = RIPDemiModel.getRIPStatic(Ref,std,TrY(i,:)',M);
               EstRIP(i) = out;
            end
            obj.EstRIP = EstRIP;
            
            % check if a sign flip occured
            if sign(obj.XTrCor)==sign(obj.XLOOCor), obj.Flipped = false; return; end
            disp('Flipping');
            obj.EstRIP = -1*obj.EstRIP;
            obj.Flipped = true;
            %toc;
            %obj.EstRIP = getRIP(obj,obj.TrY');
            %obj.EstRIP = obj.TrRIP';
        end
        function deployModel(obj)
           demiModel.deployModel(obj);
           %obj.TrRIP = [];
        end
    end
    methods (Static = true)
        function out = getRIPStatic(Ref,std,in,M)
            [~,n] = size(in); % determine input size
            out = nan*zeros(1,n);% allocate memory
            % distance between input faces and used reference
            Dir = repmat(Ref,1,n);
            Dir = (in-Dir);
            dist = sqrt(sum(((Dir./repmat(std,1,n)).^2)));
            % Direction between input face and used reference
            coeff2 = M(end,:)'./std;
            for i=1:1:n
                coeff1 = Dir(:,i)/norm(Dir(:,i));
                % Direction of RIP model
                coeff1 = coeff1./std;
                % Angle between direction and model direction
                T = coeff1'*coeff2;
                N = sqrt((coeff1'*coeff1)*(coeff2'*coeff2));
                angle = T/N;
                % parallell distance (RIP)
                out(i) = angle*dist(i);
            end   
        end
    end
end


%             estRIP = zeros(size(obj.TrRIP));
%             parfor i=1:obj.nrTr
%                 estRIP(i) = getRIP(obj,obj.TrY(i,:)');
%             end
%             obj.EstRIP = estRIP;

%             [~,n] = size(in); % determine input size
%             out = nan*zeros(1,n);% allocate memory
%             std = obj.YEigStd;
%             % distance between input faces and used reference
%             Dir = repmat(obj.YAvg,1,n);
%             Dir = (in-Dir);
%             dist = sqrt(sum(((Dir./repmat(std,1,n)).^2)));
%             % Direction between input face and used reference
%             coeff2 = obj.TrM(end,:)'./std;
%             for i=1:1:n
%                 coeff1 = Dir(:,i)/norm(Dir(:,i));
%                 % Direction of RIP model
%                 coeff1 = coeff1./std;
%                 % Angle between direction and model direction
%                 T = coeff1'*coeff2;
%                 N = sqrt((coeff1'*coeff1)*(coeff2'*coeff2));
%                 angle = T/N;
%                 % parallell distance (RIP)
%                 out(i) = angle*dist(i);
%             end
