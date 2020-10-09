classdef EvoMorphMV < superClass
    % CONTAINER PROPERTIES
    properties (Transient = true)% To make sure that certain object are not co-deleted, co-cloned or co-saved
       PPMData;% Container containing all the SNP information
       ModData;% Container containing Covariates and Shape Space Informtation
       X = [];
       FaceSpace = [];
       EvalInd;
       LCL = [];
    end
    properties (Hidden = true) 
       NormFaceSpace = false;% Normalized (mahalanobis) ShapeSpace or not
       Adjust = 'RESCALE'; % in case you perform optimization in normalized space
       Trim = 2;
       MovieFile = [];
       MovieRecord = false;
    end
    properties (Hidden = true, Dependent = true)
       RefFace;            % Reference Face
       AvgC;               % Reference Coefficients
       StdC;               % Eigenvalues of Face Space
       StdO;               % Eigenvalues to use during optimization
       nrCov;              % Number of Covariates
       nrGB;               % Number of Genetic Background Axes
       nrBF;               % Number of Covariates + Genetic Background Axes
       UseInd;
    end
    % ALGORITHMIC PROPERTIES
    properties 
       Pop; % Current Population
       PopSize = 2000; % size of the Population
       GenerationGenerator = 'NP'; % Method to generate next generation NORM (Parametric Normal Distribution), NP (Non Parametric Distribution), CR (Crossover and Mutation, traditional genetic algorithm
       Reduction = 0.5; % Fraction to reduce
       ScalingFcn = 'Rank';% Scaling options
       ResFaceType = 'WeightedAverage';% WeightedAverage Average Best
       ResFace = [];% resulting face
       ResFaceC = [];% resulting face normalized coefficients
       WFace = [];
       WFaceC = [];
       WFaceAtyp = [];
       BFace = [];
       BFaceC = [];
       BFaceAtyp = [];
       WFaceM = [];
       WFaceCM = [];
       WFaceAtypM = [];
       BFaceM = [];
       BFaceCM = [];
       BFaceAtypM = [];
       ResFaceSpace = [];% If multiple evolution are performed, the output is a space stored in this property
       Scores = 0;% scores for all Pop individuals
       InitFcn = 'Normal';% uniform in ranges, Gaussian according to eigenvalues,Training given populatio
       MaxGenerations = 100;
       Display = false;
       Verbose = false;
       EvolveCycles = 'Single';% the number of evolutions to perform: single or multiple
       nrCycles = 100;% Number of cycles to perform in case Evolve Cycles is set to multiple 
       EvalMatch = 'DET'; % Stoch, Stochastic; Det, Deterministic
       nrStoch = 100;% number random each run;
       FaceReg = false;% FACE Regularization on or off
       FKappa = 3;% Determines the range of plausible facial solutions
    end
    properties (Dependent = true)
       Diversity;pFaces;
       BestScore;MeanScore;WorstScore;
       PopAverage;
       RedInd;
       FL;
    end
    properties (Hidden = true)
       GlobalScale = 3;% Global search space factor times standard deviations in face model
       LocalScale = 1;% Local search space factor times standard deviations in face model 
       State = 'Init';
       MaxTime = +inf;
       FLimit = -inf;
       StallGenLimit = 50;
       FunctionTolerance = 0.0001;   
       TrimPop = false;
       Generation = 0;
       StallGeneration = 0;
       Time = 0;
       StartTime = 0;
       GenoTypes = [];
       scaledScores = [];
       OldBestScore = +inf;
       OldPop = [];
       OldGeneration = nan;
       OldScores = [];
    end
    properties (Hidden = true, Dependent = true)
       Dim;
       BestInd;
       WorstInd;
       sortIndex;
       GlobalRange;
       LocalRange;
    end
    properties (Hidden = true)
       Tol = 10;
       mu = 0;
       Sigma = 0;
    end
    % VISUALISATION
    properties (Hidden = true)
       Fig = [];
       sError = [];
       sDiversity = [];
       sScores = [];
       sScaling = [];
       vFace = [];
       vPop = [];
       PopLM = [];
       vAxes = [1 2 3];
       nrNoImprove = 0;
    end
    methods % CONSTRUCTOR
        function obj = EvoMorphMV(varargin)
            obj = obj@superClass(varargin{:});
        end
    end
    methods % GENERAL GETTING/SETTING
        function out = get.PPMData(obj)
           out = obj.PPMData;
           if ~superClass.isH(out), out = []; end
        end
        function out = get.Diversity(obj)
            out = mean(sqrt(sum((obj.Pop-repmat(obj.PopAverage,obj.PopSize,1)).^2)));
        end
        function out = get.BestScore(obj)
            out = min(obj.Scores);
        end
        function out = get.MeanScore(obj)
            out = mean(obj.Scores);
        end
        function out = get.WorstScore(obj)
            out = max(obj.Scores);
        end
        function out = get.BestInd(obj)
            [~,out] = min(obj.Scores);
        end
        function out = get.WorstInd(obj)
            [~,out] = max(obj.Scores);
        end
        function out = get.PopAverage(obj)
            out = mean(obj.Pop);
        end
        function out = get.Dim(obj)
            if isempty(obj.FaceSpace), out = 0; return; end
            out = obj.FaceSpace.nrEV;
        end
        function out = get.AvgC(obj)
           if isempty(obj.FaceSpace), out = 0; return; end
           out = obj.FaceSpace.AvgCoeff;
        end
        function out = get.StdC(obj)
           if isempty(obj.FaceSpace), out = 0; return; end
           out = obj.FaceSpace.EigStd;
        end  
        function out = get.GlobalRange(obj)
           out = [obj.AvgC-obj.GlobalScale*obj.StdO obj.AvgC+obj.GlobalScale*obj.StdO];
        end
        function out = get.LocalRange(obj)
            out = [obj.AvgC-obj.LocalScale*obj.StdO obj.AvgC+obj.LocalScale*obj.StdO];
        end
        function out = get.sortIndex(obj)
            if isempty(obj.Scores), out = []; return; end
            [~,out] = sort(obj.Scores);
        end
        function out = get.PopLM(obj)
            out = obj.PopLM;
            if ~superClass.isH(out), out = []; end
        end
        function out = get.pFaces(obj)
            if isempty(obj.FaceSpace), out = []; return; end
            out = sum(((obj.Pop-repmat(obj.AvgC',obj.PopSize,1))./repmat(obj.StdO',obj.PopSize,1)).^2,2)'/obj.Dim;
        end
        function out = get.RefFace(obj)
                 if isempty(obj.FaceSpace), out = []; return; end
                 out = obj.FaceSpace.Average;
        end
        function out = get.StdO(obj)
                 if obj.NormFaceSpace, out = ones(obj.Dim,1); return; end
                 out = obj.StdC;
        end
        function out = get.RedInd(obj)
                 out = ceil(obj.Reduction*obj.PopSize); 
        end
        function out = get.FL(obj)
           %out = (1./sqrt(((2*pi)^2).*obj.StdO)).*exp(-0.5*obj.FKappa^2);
           out = (1./(sqrt(2*pi).*obj.StdO)).*exp(-0.5*obj.FKappa^2);
        end
        function out = get.UseInd(obj)
           if isempty(obj.X), out = []; return; end
           tmp = obj.X(obj.EvalInd);
           out = obj.EvalInd(~isnan(tmp));
        end
    end
    methods % MAIN INTERFACE FUNCTIONS
        function evolveSingle(obj)
            obj.State = 'Init';
            resetFaces(obj);
            if obj.MovieRecord&&obj.Display
               obj.MovieFile = VideoWriter('evomovie.avi','Uncompressed AVI');
               open(obj.MovieFile);
            end
            while ~strcmp(obj.State,'Stop');
                if strcmp(obj.State,'Init')
                    obj.StartTime = cputime; obj.Time = 0;obj.Generation = 1;
                    InitializePop(obj);
                else
                    obj.Generation = obj.Generation + 1;
                    obj.Time = cputime-obj.StartTime;
                    nextGeneration(obj);
                end
                evalPop(obj);
                scaleScores(obj); 
                if obj.Display 
                   plotEvolution(obj);
                   if obj.MovieRecord
                      f = getframe(obj.vFace.Figure);
                      writeVideo(obj.MovieFile,f);
                   end
                end
                stopTest(obj);
                %pause;
            end
            obj.State = 'Final';
            if obj.MovieRecord&&obj.Display
               close(obj.MovieFile);
            end
            %if obj.SolveLocal,solveLocally(obj);end
            %beep;
        end
        function InitializePop(obj)
            if isempty(obj.FaceSpace), error('First Load Face Model'); end
            switch obj.InitFcn
                case 'Uniform'
                    R = obj.GlobalRange;
                    obj.Pop = (repmat(R(:,1),1,obj.PopSize) + repmat(R(:,2)-R(:,1),1,obj.PopSize).*rand(obj.Dim,obj.PopSize))';
                case 'Normal'
                    obj.Pop = (repmat(obj.AvgC,1,obj.PopSize) + repmat(obj.StdO,1,obj.PopSize).*randn(obj.Dim,obj.PopSize))';
                case 'Training'
                    ind = randsample(obj.FaceSpace.n,obj.PopSize,true);
                    obj.Pop = obj.FaceSpace.Tcoeff(ind,:);
                    if obj.NormFaceSpace, obj.Pop = obj.Pop./repmat(obj.StdC',length(ind),1); end
                otherwise
                    error('WRONG INITIALIZATION PROCEDURE');
            end
            trimPop(obj);
        end
        function evalPop(obj)
            switch obj.EvalMatch
                case 'DET'
                    if strcmp(obj.State,'Init')
                       obj.Scores = evalSolutions(obj,obj.Pop./repmat(obj.StdO',obj.PopSize,1));
                    else
                       red = obj.RedInd; 
                       obj.Scores = obj.Scores(obj.sortIndex);
                       redPop = obj.Pop(red+1:end,:);
                       obj.Scores(red+1:end) = evalSolutions(obj,redPop./repmat(obj.StdO',size(redPop,1),1));
                    end
                case 'STOCH'
                    if strcmp(obj.State,'Final')
                       obj.Scores = evalSolutions(obj,obj.Pop./repmat(obj.StdO',obj.PopSize,1));
                    else
                       obj.Scores = evalSolutionsStoch2(obj,obj.Pop./repmat(obj.StdO',obj.PopSize,1));
                    end
            end
        end        
        function nextGeneration(obj)
            switch obj.GenerationGenerator
                case 'NORM'
                    NORMNextGeneration(obj);
                case 'NP'
                    NPNextGeneration(obj);
                case 'CR'
                    CRNextGeneration(obj);
            end
            trimPop(obj);
        end
        function solveLocally(obj)
            switch obj.ResFaceType
                case 'WeightedAverage'
                    X0 = sum(repmat(obj.scaledScores,1,obj.Dim).*obj.Pop)/sum(obj.scaledScores);
                case 'Average'
                    X0 = obj.PopAverage;
                case 'Best'
                    X0 = obj.Pop(obj.BestInd,:);
            end
            options = optimset('Display','Iter');
            switch obj.LocalSolver
                case 'fminsearch'
                    [obj.LocalSolution,obj.LocalScore] = fminsearch(@(X) myLocalSolverFun(obj,X),X0,options);
                case 'patternsearch'
                    [obj.LocalSolution,obj.LocalScore] = patternsearch(@(X) myLocalSolverFun(obj,X),X0,[],[],[],[],[],[],[],options);
            end
        end
    end
    methods % AID INTERFACE FUNCTIONS
        function linkWithData(obj,COV,GB,SNP,RS)
            obj.X = getTestX(obj.PPMData,COV,GB,SNP,RS);
        end
        function NORMNextGeneration(obj)
                 red = floor(obj.Reduction*obj.PopSize);
                 obj.mu = mean(obj.Pop(obj.sortIndex(1:red),:));
                 obj.Sigma = cov(obj.Pop(obj.sortIndex(1:red),:));
                 obj.Pop = mvnrnd(obj.mu,obj.Sigma,obj.PopSize);
        end
        function NPNextGeneration(obj)
                 red = obj.RedInd;
                 obj.Pop = obj.Pop(obj.sortIndex,:);
                 PopR = obj.Pop(1:red,:);
                 RemainSize = obj.PopSize - red;
                 h = 1.06*min(iqr(PopR)/1.349,std(PopR))./(red^(1/5));
                 PopR = PopR';c= size(PopR,2);
                 randsel = bsxfun(@(x,y) x(randperm(y(1))),PopR',repmat(c,1,c)');
                 obj.Pop(red+1:end,:) = randsel(1:RemainSize,:) + repmat(h,RemainSize,1).*randn(RemainSize,obj.Dim); 
                 %obj.Pop(red+1:end,:) = bsxfun(@(x,y) x(randperm(y(1))),PopR',repmat(c,1,c)') + repmat(h,RemainSize,1).*randn(RemainSize,obj.Dim); 
        end
        function CRNextGeneration(obj)
            % OLD SCHOOL GENETIC ALGORITM
            % TO BE IMPLEMENTED
            obj.RedInd;
            return;
        end
        function scaleScores(obj)
            switch obj.ScalingFcn
                case 'Rank'
                    expectation = fitscalingrank(obj);
                case 'Top'
                    expectation = fitscalingtop(obj);
                case 'Shift'
                    expectation = fitscalingshiftlinear(obj);
                case 'Proportional'
                    expectation = fitscalingprop(obj);
            end
            obj.scaledScores = expectation;
        end
        function stopTest(obj)
            obj.State = 'Run';
            if obj.Generation>obj.MaxGenerations, obj.State = 'Stop'; end
            if obj.Time>obj.MaxTime, obj.State = 'Stop'; end
        end
        function trimPop(obj)
            % Pop must remain within global range
            if ~obj.TrimPop, return; end
            M = repmat(obj.GlobalRange(:,1)',obj.PopSize,1);
            index = find(obj.Pop<M);
            if ~isempty(index),obj.Pop(index) = M(index);end
            M = repmat(obj.GlobalRange(:,2)',obj.PopSize,1);
            index = find(obj.Pop>M);
            if ~isempty(index),obj.Pop(index) = M(index);end
        end
        function out = getWAverage(obj)
            W = obj.scaledScores;
            W(obj.sortIndex(obj.RedInd:end)) = 0;
            out = sum(repmat(W,1,obj.Dim).*obj.Pop)/sum(W);
            %out = sum(repmat(obj.scaledScores,1,obj.Dim).*obj.Pop)/sum(obj.scaledScores);
        end
        function resetFaces(obj)
            obj.ResFace = clone(obj.RefFace);
            obj.WFace = [];obj.WFaceC = [];obj.WFaceAtyp = [];
            obj.BFace = [];obj.BFaceC = [];obj.BFaceAtyp = [];
        end
    end
    methods % MATCHING FUNCTIONS
        function out = evalSolutions(obj,sol)
            DepVar = projectSolutionsFast(obj,sol);
            % MATCHING
            out = matchSolutions(obj.PPMData,obj.X,DepVar,obj.UseInd);
            if ~obj.FaceReg, return; end
            % PROBABILITY FACES
            FM = mean(matchF(obj,sol),1);
            out = (out+FM)/2;
            % Facial Plausibilities
            %Matches = [Matches;matchF(obj,sol)];
            % COMBINING
            %out = nansum(repmat(obj.MW,1,n).*Matches,1)./repmat(nansum(obj.MW),1,n);
        end
        function out = matchF(obj,sol,ind)
            if nargin < 3, ind = (1:obj.Dim);end
            n = size(sol,1);
            out = normpdf(sol',repmat(obj.AvgC,1,n),repmat(obj.StdO,1,n));
            out = -log(out./(out+repmat(obj.FL,1,n)));
            out = out(ind,:);
        end
        function out = projectSolutions(obj,sol)
             n = size(sol,1);
            % RECONSTRUCT FACES
            faces = repmat(obj.FaceSpace.AvgVec,1,n) + obj.FaceSpace.EigVec*sol';
            faces = reshape(faces,3,obj.FaceSpace.Average.nrV,n);
            % PROJECT FACES
            LCL = obj.LCL;
            Input = cell(1,length(LCL));
            SD = obj.ModData;
            f = waitbar(0,'PROJECTIONS');
            counter = 0;
            for l=1:1:length(LCL)
                % l=1;
                Input{l}.DepVar = cell(1,LCL(l));
                for i=1:1:LCL(l)
                    counter = counter+1;
                    % i=1;
                    space = SD{l}.SymSpace{i};
                    spaceind = SD{l}.SubInd{i};
                    tmp = nan*zeros(n,space.nrEV);
                    % crop
                    forfaces = faces(:,spaceind,:);
                    parfor j=1:1:n
                        %j = 1;
                        % align
                        T = scaledRigidTM;
                        points = squeeze(forfaces(:,:,j));
                        match(T,space.AvgVertices,points);
                        eval(T,points); %#ok<PFBFN>
                        D = T.Evaluation.Vertices(:);
                        % project
                        D = Vec2Coeff(space,D);
                        D = D./space.EigStd;
                        tmp(j,:) = D;
                    end
                    Input{l}.DepVar{i} = tmp;
                    waitbar(counter/sum(LCL),f);
                end
            end
            delete(f);
            out = Input;
        end
        function out = projectSolutionsFast(obj,sol)
            n = size(sol,1);
            % RECONSTRUCT FACES
            faces = repmat(obj.FaceSpace.AvgVec,1,n) + obj.FaceSpace.EigVec*sol';
            faces = reshape(faces,3,obj.FaceSpace.Average.nrV,n);
            % PROJECT FACES
            LCL = obj.LCL;
            Input = cell(1,length(LCL));
            SD = obj.ModData;
            f = waitbar(0,'PROJECTIONS');
            counter = 0;
            for l=1:1:length(LCL)
                % l=1;
                Input{l}.DepVar = cell(1,LCL(l));
                for i=1:1:LCL(l)
                    counter = counter+1;
                    % i=1;
                    space = SD{l}.SymSpace{i};
                    spaceind = SD{l}.SubInd{i};
                    % crop
                    forfaces = faces(:,spaceind,:);
                    forref = space.AvgVertices;
                    newforfaces = nan*zeros(3*size(forref,2),n);
                    % superimpose
                    parfor j=1:n
                        %[~,Z] = procrustes(forref',squeeze(forfaces(:,:,j))','reflection',false);
                        points = squeeze(forfaces(:,:,j));
                        T = EvoMorphMV.posematch(forref,points);
                        points = T.Scale*(T.Rotation*(points)+repmat(T.Translation,1,size(points,2)));
                        newforfaces(:,j) = points(:);
                    end
                    Input{l}.DepVar{i} = (space.EigVec'*(newforfaces-repmat(space.AvgVec,1,n)))';
                    waitbar(counter/sum(LCL),f);
                end
            end
            delete(f);
            out = Input;
        end
    end
    methods % Plotting functions
        function updateResFace(obj)
            switch obj.ResFaceType
                case 'WeightedAverage'
                    [out,C] = getWAverageFace(obj);
                case 'Average'
                    [out,C] = getAverageFace(obj);
                case 'Best'
                    [out,C] = getBestFace(obj);
            end
            obj.ResFace.Vertices = out.Vertices;
            obj.ResFaceC = C;
        end
        function updatePopLM(obj)
            obj.PopLM.Vertices = obj.Pop(:,obj.vAxes)';
            obj.PopLM.Value = obj.Scores;
        end
        function plotEvolution(obj)
            if strcmp(obj.State,'Init');
                constructDisplay(obj);
            end
            updateResFace(obj);
            updatePopLM(obj);
            plot(obj.sError,obj.Generation,obj.BestScore,'b.','MarkerSize',10);
            plot(obj.sError,obj.Generation,obj.MeanScore,'k.','MarkerSize',10);
            title(obj.sError,['Best score: ' num2str(obj.BestScore) ' Mean Score: ' num2str(obj.MeanScore)]);
            %title(obj.sError,[' Error Score: ' num2str(obj.MeanScore)]);
            plot(obj.sDiversity,obj.Generation,obj.Diversity,'k.','MarkerSize',10);
            hist(obj.sScores,obj.Scores,20);
            set(obj.sScores,'xlim',[obj.BestScore obj.WorstScore],'ylim',[0 obj.PopSize/2]);
            plot(obj.sScaling,obj.Scores,obj.scaledScores,'b.');
            drawnow;
        end
    end
    methods % FACIAL SOLUTIONS
        function [out,C,Atyp] = getWorstFace(obj)
            C = obj.Pop(obj.WorstInd,:);
            [out,Atyp] = getFaceFromC(obj,C);
        end
        function [out,C,Atyp] = getAverageFace(obj)
            C = obj.PopAverage;
            [out,Atyp] = getFaceFromC(obj,C);
        end
        function [out,C,Atyp] = getBestFace(obj)
            C = obj.Pop(obj.BestInd,:);
            [out,Atyp] = getFaceFromC(obj,C);
        end
        function out = getLocalFace(obj)
            out = getScan(obj.FaceModel,obj.LocalSolution);
            %out = getFaceFromC(obj,C);
        end
        function [out,C,Atyp] = getWAverageFace(obj)
            C = getWAverage(obj);
            [out,Atyp] = getFaceFromC(obj,C);
        end
        function [out,Atyp] = getFaceFromC(obj,C)
                 if ~obj.NormFaceSpace, out = getScan(obj.FaceSpace,C); Atyp = sqrt(sum((C./obj.StdC').^2)); return; end
                 Atyp = (sum(C.^2))/length(C);
                 switch obj.Adjust
                     case 'RESCALE'
                         C = C.*obj.StdC';
                     case 'TRIM'
                         B = obj.SOTrim*obj.StdC';
                         index = find(abs(C)>B);
                         C(index) = B(index).*sign(C(index));
                     otherwise
                 end
                 out = getScan(obj.FaceSpace,C);
        end
        function out = getSolutionSpace(obj)
            out = clone(obj.FaceModel);
            out.Tcoeff = obj.Pop;
            Data = reconstructTraining(out);
            switch out.Type
                case 'shapePCA'
                    getAverage(out,Data);
                    getModel(out,Data);
                case 'appearancePCA'
                    shapePC = Data(1:out.nrSC,:)/out.WS;
                    out.Shape.Tcoeff = shapePC';
                    dataShape = reconstructTraining(out.Shape);
                    outShape = clone(out.Shape);
                    getAverage(outShape,dataShape);
                    getModel(outShape,dataShape);
                    clear dataShape shapePC;
                    texPC = Data(out.nrSC+1:end,:);
                    out.Texture.Tcoeff = texPC';
                    dataTex = reconstructTraining(out.Texture);
                    outTex = clone(out.Texture);
                    getAverage(outTex,dataTex);
                    getModel(outTex,dataTex);
                    out.Shape = outShape;
                    out.Texture = outTex;
                    getModel(out);
            end
            
            obj.ResFaceSpace = out;
        end
    end
    methods % Other general interface functions
        function delete(obj)
            try
                close(obj.Fig);
            catch
            end
            %superClass.delete(obj);
        end
    end
    methods % Main Prediction function + amplifying and exporting functions
        function [scan,Atyp,refscan] = getFacialPrediction(obj)
                 %initializeInput(obj);
                 evolve(obj);
                 updateResFace(obj);
                 %[scan,~,Atyp] = getBestFace(obj);
                 [scan,~,Atyp] = getWAverageFace(obj);
                 %inVec = getVec(obj.FaceModel,obj.Subject.BaseFace);  
                 %C = weightedFit2(obj.FaceModel,inVec,[],0.5,(1:length(inVec)));
                 %C = getCoeff(obj.FaceModel,obj.Subject.BaseFace);
                 %refscan = getScan(obj.FaceModel,C);
                 refscan = obj.RefFace;
        end
    end
    methods (Static = true)
       function out = posematch(q,p)
        % find the rigid transform from p to q
                 nbpts = size(p,2);
                 % align centers
                 centerq = mean(q,2);
                 centerp = mean(p,2);
                 q = q-repmat(centerq,1,nbpts);
                 p = p-repmat(centerp,1,nbpts);
                 % Scale both by making mean length 1
                 lengths = sqrt(sum(p.^2,1));
                 meanLength1 = mean(lengths,2);
                 p = p./meanLength1;
                 lengths = sqrt(sum(q.^2,1));
                 meanLength2 = mean(lengths,2);
                 q = q./meanLength2;
                 scalefactor = meanLength2/meanLength1;
                 % Rotate
                 [U,S,V] = svd(q*p');
                 H = V*sign(S)*U';
                 H = H';
                 % putting it all together
                 out.Scale=scalefactor;
                 out.Rotation = H;
                 transform = out.Scale*H;
                 Ta=eye(4);
                 Ta(1:3,4)=centerq;%+obj.c;
                 Tb=eye(4);
                 Tb(1:3,4)=-centerp;%-obj.c;
                 R=eye(4);
                 R(1:3,1:3)=transform;
                 Tout=Ta*R*Tb;
                 out.Translation = Tout(1:3,4)/out.Scale;
       end
    end
end


%for i=1:obj.Dim
                 %   [~,~,h] = ksdensity(PopR(:,i));
                 %   mark = randsample(red,RemainSize,true);
                 %   obj.Pop(red+1:end,i) = PopR(mark,i)+h.*randn([RemainSize 1]);                
                 %end

%if obj.NormShapeSpace, out = sum((obj.Pop-repmat(obj.AvgC',obj.PopSize,1)).^2,2)/obj.Dim; return; end

%if isempty(obj.FaceSpace), out = []; return; end
%if obj.NormShapeSpace, out = [obj.AvgC-obj.LocalScale*ones(obj.Dim,1) obj.AvgC+obj.LocalScale*ones(obj.Dim,1)]; return; end
            

%if isempty(obj.FaceSpace), out = []; return; end
%if obj.NormShapeSpace, out = [obj.AvgC-obj.GlobalScale*ones(obj.Dim,1) obj.AvgC+obj.GlobalScale*ones(obj.Dim,1)]; return; end
            

%             obj.Scores = (1-obj.FacePenalty)*...
%                          matchFaces(obj.SNPCont,obj.Pop./repmat(obj.FaceModel.EigStd',obj.PopSize,1),obj.GenoTypes,obj.BF)...
%                          + obj.FacePenalty*obj.pFaces;


%         function out = get.BF(obj)
%                  if isempty(obj.Subject), out = 0; return; end
%                  out = obj.Subject.BF;
%         end

%         function out = get.Subject(obj)
%            out = obj.Subject;
%            if ~superClass.isH(out), out = []; end
%         end


%         function out = get.FaceSpace(obj)
%            if isempty(obj.BaseCont), out = []; return; end
%            if ~obj.UseRedShape, out = obj.BaseCont.ShapeSpace; return; end 
%            out = obj.BaseCont.RedShapeSpace;
%         end


%             
%             obj.MatchW
%             
%             if obj.MatchUse==0
%                 BaseMatches = -log(matchFaces(obj.BaseCont,sol,obj.CovDistr,obj.GBDistr)); 
%                 if strcmp(obj.EvalMatch,'Stoch')
%                     ind = randsample(obj.nrSNP,obj.nrR);
%                     ind = sort(ind);
%                     diffind = setdiff((1:obj.nrSNP),ind);
%                     Wind = setdiff(1:obj.nrBF+obj.nrSNP+obj.Dim,diffind+obj.nrBF);
%                 else
%                     ind = (1:obj.nrSNP);
%                     Wind = (1:obj.nrBF+obj.nrSNP+obj.Dim);
%                 end
%                 SNPMatches = -log(matchFaces(obj.SNPCont,sol,obj.GT,obj.BF,ind));
%                 Matches = [BaseMatches;SNPMatches];
%             elseif obj.MatchUse==1
%                 Matches = -log(matchFaces(obj.BaseCont,sol,obj.CovDistr,obj.GBDistr));
%                 Wind = (1:obj.nrBF+obj.Dim);
%             else
%                 Matches = -log(matchFaces(obj.SNPCont,sol,obj.GT,obj.BF));
%                 Wind = (1:obj.nrSNP+obj.Dim);
%             end
%             pMatches = -log(FMatches(obj,sol));
%             out = nansum(repmat(obj.MatchW(Wind),1,n).*[Matches;pMatches],1)./repmat(nansum(obj.MatchW(Wind)),1,n);


% function evolve(obj)
%             switch obj.EvolveCycles
%                 case 'Single'
%                     evolveSingle(obj);
%                 case 'Multiple'
%                     obj.Display = false;
%                     CyclePop = nan*zeros(obj.Dim,obj.nrCycles);
%                     warning off;
%                     tic;
%                     %parfor_progress(obj.nrCycles);
%                     FM = obj.FaceModel;
%                     DM = obj.DemiModelSet;
%                     parfor i=1:obj.nrCycles
%                         out = [];
%                         forobj = clone(obj);
%                         forobj.FaceModel = FM;% passing the pointer not the whole object to save memory (see transient aspect)
%                         forobj.DemiModelSet = DM;
%                         forobj.EvolveCycles = 'Single';
%                         forobj.Display = false;
%                         evolve(forobj);
%                         switch obj.ResFaceType
%                             case 'WeightedAverage'
%                                 out = sum(repmat(forobj.scaledScores,1,forobj.Dim).*forobj.Pop)/sum(forobj.scaledScores);
%                             case 'Average'
%                                 out = forobj.PopAverage;
%                             case 'Best'
%                                 out = forobj.Pop(forobj.BestInd,:);
%                         end
%                         CyclePop(:,i) = out(:);
%                         delete(forobj);
%                         %parfor_progress;
%                     end
%                     toc;
%                     warning on;
%                     obj.PopSize = obj.nrCycles;
%                     obj.Pop = CyclePop';
%                     evalPop(obj);
%                     scaleScores(obj);
%                     getSolutionSpace(obj);
%                     %parfor_progress(0);
%             end
%             updateResFace(obj);
%         end


%     function out = evalSolutionsStoch(obj,sol)
%             n = size(sol,1);
%             % DNA Matches
%             if ~strcmp(obj.MatchWhat,'ALL'), out = evalSolutions(obj,sol); return; end
%             ind = sort(randsample(obj.indT,obj.nrStoch));
%             Matches = [matchBF(obj,sol,intersect(obj.indBF,ind));...
%                        matchSNP(obj,sol,intersect(obj.indSNP,ind)-obj.nrBF);...
%                        matchF(obj,sol,intersect(obj.indDIM,ind)-obj.nrBF-obj.nrSNP)];
%             % COMBINING
%             out = nansum(repmat(obj.MW(ind),1,n).*Matches,1)./repmat(nansum(obj.MW(ind)),1,n);
%         end
%         function out = evalSolutionsStoch2(obj,sol)
%             n = size(sol,1);
%             % DNA Matches
%             if ~strcmp(obj.MatchWhat,'ALL'), out = evalSolutions(obj,sol); return; end
%             ind = sort(randsample(obj.indSNP,obj.nrStoch));
%             indW = [obj.indBF, ind, obj.indDIM];
%             Matches = [matchBF(obj,sol);...
%                        matchSNP(obj,sol,ind-obj.nrBF);...
%                        matchF(obj,sol)];
%             % COMBINING
%             out = nansum(repmat(obj.MW(indW),1,n).*Matches,1)./repmat(nansum(obj.MW(indW)),1,n);
%         end
%         function out = matchBF(obj,sol,ind)
%             if nargin < 3, ind = (1:obj.nrBF);end
%             switch obj.BFMatch
%                    case 'IB'
%                       out = matchFaces(obj.BaseCont,sol,obj.CovDistr,obj.GBDistr);
%                    case 'pTIB'
%                       [~,out] = matchFaces(obj.BaseCont,sol,obj.CovDistr,obj.GBDistr);
%                    case 'pT'
%                       [~,~,out] = matchFaces(obj.BaseCont,sol,obj.CovDistr,obj.GBDistr);
%             end
%             out = -log(out);
%             out = out(ind,:);
%         end
%         function out = matchSNP(obj,sol,ind)
%             if nargin < 3, ind = (1:obj.nrSNP);end
%             switch obj.SNPMatch
%                    case 'IB'
%                       out = matchFaces(obj.SNPCont,sol,obj.GT,obj.BF,ind);
%                    case 'pTIB'
%                       [~,out] = matchFaces(obj.SNPCont,sol,obj.GT,obj.BF,ind);
%                    case 'pT'
%                       [~,~,out] = matchFaces(obj.SNPCont,sol,obj.GT,obj.BF,ind);
%             end
%             out = -log(out);
%         end
%         
