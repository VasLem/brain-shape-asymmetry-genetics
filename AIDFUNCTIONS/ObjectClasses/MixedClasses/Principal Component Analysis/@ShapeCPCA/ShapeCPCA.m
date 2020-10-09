classdef ShapeCPCA < superClass
    properties
        X = []; % Data
    end
    %% Super PCA
    properties % Super PCA
        nrEV;
        Scale = false;
    end
    properties
        RefScan = [];
    end
    properties (Dependent = true)
        Scores; % Super Scores
        EigStd;% Standard deviations per component
        Explained;% Percentage Explained in Eigenvectors
        n;% number of trainings Data
        nrLandmarks;
    end
    properties ( Hidden = true )
        U = [];% U  in SVD
        S = [];% s in SVD
        V = [];% V in SVD
    end
    %% Block PCA
    properties % Block PCA
        BlockScale = false;
    end
    properties
        BlockLabel;
        nrBlocks;
        BlockSizes;
        BlockScores; % Block
        BlockLoadings;
        BlockWeights; % Block contribution to the super scores
        BlockExplained; % Explained Variance in each Block
    end
    methods % Constructor
        function obj = ShapeCPCA(varargin)
            obj = obj@superClass(varargin{:});
        end
    end
    methods % Get & Set
        function out = get.n(obj)
            if isempty(obj.X), out = []; return; end;
            out = size(obj.X,1);
        end
        function out = get.Scores(obj)
            if isempty(obj.X), out = []; return; end;
            if isempty(obj.U) || isempty (obj.S) || isempty (obj.V)
                if obj.Scale
                    Z = zscore(obj.X);
                else
                    meanX = mean(obj.X);
                    Z = bsxfun(@minus,obj.X,meanX);
                end
                [obj.U,obj.S,obj.V] = svd(Z,0);
            end
            out = obj.U(:,1:obj.nrEV)*obj.S(1:obj.nrEV,1:obj.nrEV);
        end
        function out = get.EigStd(obj)
            if isempty(obj.X), out = []; return; end;
            if isempty(obj.U) || isempty (obj.S) || isempty (obj.V)
                if obj.Scale
                    Z = zscore(obj.X);
                else
                    meanX = mean(obj.X);
                    Z = bsxfun(@minus,obj.X,meanX);
                end
                [obj.U,obj.S,obj.V] = svd(Z,0);
            end
            Std = diag(obj.S);
            out = Std(1:obj.nrEV)./sqrt(obj.n-1);
        end
        function out = get.Explained(obj)
            if isempty(obj.X), out = []; return; end;
            if isempty(obj.U) || isempty (obj.S) || isempty (obj.V)
                if obj.Scale
                    Z = zscore(obj.X);
                else
                    meanX = mean(obj.X);
                    Z = bsxfun(@minus,obj.X,meanX);
                end
                [obj.U,obj.S,obj.V] = svd(Z,0);
            end
            Var = diag(obj.S);
            Var = Var.^2./(obj.n-1);
            out = sum(obj.EigStd.^2)/sum(Var);
        end
        function out = get.nrLandmarks(obj)
            if isempty(obj.X), out = []; return; end;
            out = floor(size(obj.X,2)/3);
        end
        function updateBlockInfo(obj)
            if isempty(obj.LandmarksLabel)
                obj.BlockLabel = (1:obj.nrLandmarks);
                obj.BlockLabel = repmat(obj.BlockLabel,3,1);
                obj.BlockLabel = reshape(obj.BlockLabel,[1 3*obj.nrLandmarks]);
                obj.nrBlocks = obj.nrLandmarks;
                for i=1:obj.nrBlocks
                    obj.BlockSizes(i) = sum(obj.BlockLabel==i);
                end
                return;
            end
            obj.nrBlocks = length(unique(obj.LandmarksLabel));
            
            obj.BlockLabel = zeros(1,3*obj.nrLandmarks);
            for i=1:obj.nrLandmarks
                obj.BlockLabel((i-1)*3+1:i*3) = obj.LandmarksLabel(i);
            end
            
            
            
            obj.BlockSizes = zeros(1,obj.nrBlocks);
            for i=1:obj.nrBlocks
                obj.BlockSizes(i) = sum(obj.BlockLabel==i);
            end
        end
    end
    methods
        function updateModel(obj)
            obj.BlockScores = zeros(obj.nrBlocks,obj.nrEV,obj.n);
            obj.BlockLoadings = cell(1,obj.nrBlocks);
            for i=1:obj.nrBlocks
                obj.BlockLoadings{i} = zeros(obj.nrEV,obj.BlockSizes(i));
            end
            obj.BlockWeights = zeros(obj.nrBlocks,obj.nrEV);
            obj.BlockExplained = zeros(obj.nrBlocks,obj.nrEV);
            meanX = mean(obj.X);
            E = bsxfun(@minus,obj.X,meanX);
            Z = E;
            Var = zeros(1,obj.nrBlocks);
            for i=1:obj.nrBlocks
                ind = obj.BlockLabel==i;
                Var(i) = trace(Z(:,ind)'*Z(:,ind));
            end
            
            for i=1:obj.nrEV
                t_T = obj.Scores(:,i);
                T = zeros(obj.nrBlocks,obj.n);
                for j=1:obj.nrBlocks
                    ind = obj.BlockLabel==j;
                    Eb = E(:,ind);
                    pb = Eb'*t_T;
                    pb = pb./norm(pb);
                    tb = Eb*pb;
                    obj.BlockLoadings{j}(i,:) = pb;
                    obj.BlockScores(j,i,:) = tb;
                    T(j,:) = tb;
                    pb = Eb'*t_T/(t_T'*t_T);
                    E(:,ind) = E(:,ind) - t_T*pb';
                    obj.BlockExplained(j,i) = 1-trace(E(:,ind)'*E(:,ind))/Var(j);
                end
                Wt = T*t_T/(t_T'*t_T);
                Wt = Wt./norm(Wt);
                obj.BlockWeights(:,i) = Wt;
            end
        end
        function viewPC(obj,PC)
            Value = zeros(1,obj.nrLandmarks);
            for i=1:obj.nrBlocks
                if ~isempty(obj.LandmarksLabel)
                    Value(obj.LandmarksLabel==i) = obj.BlockExplained(i,PC);
                    if PC~=1
                        Value(obj.LandmarksLabel==i) = Value(obj.LandmarksLabel==i) - obj.BlockExplained(i,PC-1);
                    end
                else
                    Value(i) = obj.BlockExplained(i,PC);
                    if PC~=1
                        Value(i) = Value(i) - obj.BlockExplained(i,PC-1);
                    end
                end
                
            end
            obj.RefScan.Value = Value;
            obj.RefScan.ColorMode='Indexed';
            obj.RefScan.ViewMode='Solid'; obj.RefScan.Material='Metal'; obj.RefScan.LightMode='gouraud';
            v.SceneLightVisible = 1; v.BackgroundColor = [1 1 1];
            f1 = viewer(obj.RefScan);
        end
    end
end