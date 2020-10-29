classdef myPCA < superClass
    properties (Abstract = true)
        Average;% Average
    end
    properties(Abstract = true, Dependent = true)
        AvgVec;% Average vector description, dependent on the average
    end
    properties
        EigVal = [];% Eigen Values
        EigVec = [];% Eigen Vectors
        Tcoeff = [];% myPCA coeff of training data
        AvgCoeff = [];% The myPCA coeff of the average entity
        Centering = true;% Center data around average or not;
    end
    properties (Dependent = true)
        nrEV;% Number of Principal Components
        EigStd;% Standard deviations per component
        Explained;% Percentage Explained in Eigenvectors
        n;% number of trainings Data
        U;% U  in SVD
        S;% s in SVD
        V;% V in SVD
    end
    methods % Constructor
        function obj = myPCA(varargin)
            obj = obj@superClass(varargin{:});         
        end
    end
    methods % Special Setting & Getting
        function out = get.nrEV(obj)
            %out = length(obj.EigVal);
            out = size(obj.EigVec,2);
        end
        function out = get.EigStd(obj)
            if isempty(obj.EigVal), out = []; return; end
            out = sqrt(obj.EigVal);
        end
        function obj = set.EigStd(obj,in)
            obj.EigVal = in.^2;
        end
        function out = get.n(obj)
            if isempty(obj.Tcoeff), out = 0; return; end
            out = size(obj.Tcoeff,1);
        end
        function out = get.Explained(obj)
            if isempty(obj.EigVal), out = []; return; end
            out = 100*obj.EigVal/sum(obj.EigVal);
        end
        function out = get.V(obj)
            out = obj.EigVec;
        end
        function obj = set.V(obj,in)
            obj.EigVec = in;
        end
        function out = get.S(obj)
            out = sqrt(obj.EigVal)*sqrt(obj.n-1);
        end
        function obj = set.S(obj,in)
            obj.EigVal = (in./sqrt(obj.n-1)).^2;
        end
        function out = get.U(obj)
            out = obj.Tcoeff./repmat(obj.S',obj.n,1);
        end
        function obj = set.U(obj,in)
            obj.Tcoeff = repmat(obj.S',size(in,1),1).*in;
        end
        function out = get.AvgCoeff(obj)
            out = zeros(obj.nrEV,1);
        end
    end
    methods (Abstract = true) % Abstract Basic Interface & Conversion functions
        % Need to be made explicit depending on the type of data, eg shape, texture,...
        out = Struc2Vec(obj,in);% Converting Stucture representation 2 Vector representation
        out = Vec2Struc(obj,in);% Converting Vector representation 2 Structure representation
        out = IndexStruc2Vec(obj,in);% Converting Index defined on Stucture 2 Index defined on Vector representation
        %out = IndexVec2Struc(obj,in);% Converting Index defined on Vector 2 Index defined on Structure representation
        out = WeightStruc2Vec(obj,in);% Converting Weights defined on Stucture 2 Weights defined on Vector representation
        %out = WeightVec2Struc(obj,in);% Converting Weights defined on Vector 2 Weights defined on Structure representation
        out = getData(obj,in);% Simple function to extract data, from either a matrix or structure, This Data is used for model building
    end
    methods % Concrete Basic Interface & Conversion functions
       function out = Coeff2Vec(obj,in)
            % Convert given myPCA coeff into a Vector
            if size(in,1)==1, in = in'; end
            if obj.Centering
                out = obj.AvgVec + obj.EigVec*in;
            else
                out = obj.EigVec*in;
            end
       end
       function out = Vec2Coeff(obj,in,index)
            % Convert given Vector into myPCA coeff
            if nargin < 3, index = (1:length(in));end
            if length(in)~=length(index),in = in(index);end
            %length(index)
            if obj.Centering
                out = obj.EigVec(index,:)'*(in-obj.AvgVec(index));
            else
                out = obj.EigVec(index,:)'*in;
            end
        end
       function out = Coeff2Struc(obj,in)
            % Convert given myPCA coeff into a structure
            out = Vec2Struc(obj,Coeff2Vec(obj,in));
       end
       function out = Struc2Coeff(obj,in)
            % Convert a given structure into myPCA coeff
            out = Vec2Coeff(obj,Struc2Vec(obj,in));
       end
       function out = isVec(obj,in)
            % test whether input is a Vector
            out = false;
            if superClass.isH(in), return; end% if it is an object return
            [n,m] = size(in);
            if ~(n==1||m==1), return; end% if does not have a single column or single row return
            if (n==obj.nrEV||m==obj.nrEV), return; end% it is a coeff vector return
            out = true;
        end
       function out = isCoeff(obj,in)
            % test whether input is a vector of coefficients
            out = false;
            if superClass.isH(in), return; end% if it is an object return
            [n,m] = size(in);
            if ~(n==1||m==1), return; end% if does not have a single column or single row return
            if ~(n==obj.nrEV||m==obj.nrEV), return; end% it is a normal vector return
            out = true;
        end
       function out = isStruc(obj,in)
            % test whether input is a structure
            out = false;
            if isVec(obj,in)||isCoeff(obj,in), return; end
            out = true;
       end
       function out = getVec(obj,in)
            % convert input to Vector
            if isVec(obj,in),out = in;return;end
            if isCoeff(obj,in),out = Coeff2Vec(obj,in);return;end
            if isStruc(obj,in),out = Struc2Vec(obj,in);return;end 
        end
       function out = getStruc(obj,in)
            % convert input to Structure
            if isStruc(obj,in),out = in; return;end
            if isCoeff(obj,in),out = Coeff2Struc(obj,in); return;end
            if isVec(obj,in),out = Vec2Struc(obj,in);return; end
       end
       function out = getCoeff(obj,in)
            % convert input to Coeff
            if isCoeff(obj,in),out = in;return;end
            if isStruc(obj,in),out = Struc2Coeff(obj,in);return;end
            if isVec(obj,in),out = Vec2Coeff(obj,in);return;end 
       end
    end
    methods % General Interface Functions        
       function out = reconstruct(obj,in)
            % Function to create a sample based on model coeff
            out = Coeff2Struc(obj,in);
       end
       function out = fit(obj,in)
                out = getCoeff(obj,in);
       end
       function out = partialFit(obj,in,index)
            % Function to perform a basic model fit
            inVec = getVec(obj,in);
            if nargin < 3
               index = find(~isnan(inVec));
            elseif isStruc(obj,in)
               index = IndexStruc2Vec(obj,index);
            else
               %index = index;
            end
            % extract the subspace according to index
            %obj = subSpace(obj,index);
            %if length(in)~=length(index),in = in(index);end
            out = Vec2Coeff(obj,inVec,index);
            %delete(obj);
       end
       function out = weightedFit(obj,in,A,sigma,index)
             inVec = getVec(obj,in);             
             if isempty(A)||nargin<3
                A = ones(size(inVec));
             elseif isStruc(obj,in)
                A = WeightStruc2Vec(obj,A);
             end
             if isempty(sigma)||nargin<4
                sigma = 1;
             end
             if nargin < 5
                index = find(~isnan(inVec));
             elseif isStruc(obj,in)
                index = IndexStruc2Vec(obj,index);
             else
               %index = index;
             end
             if length(in)~=length(index),in = in(index);end
             if length(A)~=length(index),A = A(index);end
             A = spdiags(A,0,speye(length(in),length(in)));
             in = in-obj.AvgVec(index);
             AQ = A*obj.EigVec(index,:)*diag(obj.EigVal);
             [U,W,V] = svd(AQ,'econ');
             W = diag(W);
             W = W./((W.^2)+ones(size(W))*sigma);
             out = diag(obj.EigVal)*V*diag(W)*U'*A*in;
       end
       function out = weightedFit2(obj,in,A,sigma,index)
             inVec = in;             
             if isempty(A)||nargin<3
                A = ones(size(inVec));
             elseif isStruc(obj,in)
                A = WeightStruc2Vec(obj,A);
             end
             if isempty(sigma)||nargin<4
                sigma = 1;
             end
             if nargin < 5
                index = find(~isnan(inVec));
             elseif isStruc(obj,in)
                index = IndexStruc2Vec(obj,index);
             else
               %index = index;
             end
             if length(in)~=length(index),in = in(index);end
             if length(A)~=length(index),A = A(index);end
             A = spdiags(A,0,speye(length(in),length(in)));
             in = in-obj.AvgVec(index);
             AQ = A*obj.EigVec(index,:)*diag(obj.EigVal);
             [U,W,V] = svd(AQ,'econ');
             W = diag(W);
             W = W./((W.^2)+ones(size(W))*sigma);
             out = diag(obj.EigVal)*V*diag(W)*U'*A*in;
       end
       function out = subSpace(obj,index)
            % Function to retrieve the subspace from a myPCA model
            % subspace; a space having a subset of the variables
            if nargout == 1,obj = clone(obj);out = obj;end
            obj.AvgVec = obj.AvgVec(index);
            obj.EigVec = obj.EigVec(index,:);
       end
       function out = reduceNrPc(obj,nr)
            % Function to reduce the number of principal components
            if nargout ==1,obj = clone(obj);out = obj;end
            if nr > obj.nrEV, return; end
            obj.EigVal = obj.EigVal(1:nr);
            obj.EigVec = obj.EigVec(:,1:nr);
            obj.Tcoeff = obj.Tcoeff(:,1:nr);
       end
       function out = nrPcPercVar(obj,percvar)
            % function to determine the required number of PC for a given
            % percentage variance
            if isempty(obj.Explained),out = 0; return; end
            out = 1;perc = obj.Explained(1);
            while perc < percvar
                 perc = perc+obj.Explained(out+1);
                 out = out+1;
            end
       end
       function out = stripPercVar(obj,percvar)
           % Function to reduce the number of PC according to a given
           % percentage of variance
           if nargin < 2, percvar = 99; end
           if nargout == 1,obj = clone(obj);out = obj;end
           reduceNrPc(obj,nrPcPercVar(obj,percvar));
       end 
       function out = reconstructTraining(obj)
                out = repmat(obj.AvgVec,1,obj.n) + obj.EigVec*obj.Tcoeff';
       end
    end

end