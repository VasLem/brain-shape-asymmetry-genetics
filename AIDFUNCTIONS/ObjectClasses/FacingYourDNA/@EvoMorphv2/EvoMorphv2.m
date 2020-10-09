classdef EvoMorphv2 < superClass
    % CONTAINER PROPERTIES
    properties (Transient = true)% To make sure that certain object are not co-deleted, co-cloned or co-saved
       SNPCont;% Container containing all the SNP information
       BaseCont;% Container containing Covariates and Shape Space Informtation
    end
    properties (Hidden = true) 
       NormFaceSpace = true;% Normalized (mahalanobis) ShapeSpace or not
       Adjust = 'RESCALE'; % in case you perform optimization in normalized space
       Trim = 2;
       MovieFile = [];
       MovieRecord = false;
    end
    properties (Transient = true, Hidden = true)
       UseRedShape = false;% Kind of ShapeSpace to use
       FaceSpace = [];     % Shape Space
       GT = [];            % Genotypes
       GTFreq = [];        % Freq of Genotypes
       GTAtyp = [];        % Atypicality Genotypes
       Cov = [];           % Covariates
       GB = [];            % Genetic Background
       BF = 0;             % BF
       CovRIP = [];        % RIPPED covariates
       CovDistr = {};      % Distributions RIPPED covariates
       GBRIP = [];         % RIPPED genetic background
       GBDistr = {};       % Distributions RIPPED genetic background
       RW = [];            % Relative Weights for matching scores
    end
    properties (Hidden = true, Dependent = true)
       RefFace;            % Reference Face
       AvgC;               % Reference Coefficients
       StdC;               % Eigenvalues of Face Space
       StdO;               % Eigenvalues to use during optimization
       nrCov;              % Number of Covariates
       nrGB;               % Number of Genetic Background Axes
       nrBF;               % Number of Covariates + Genetic Background Axes
    end
    % ALGORITHMIC PROPERTIES
    properties 
       Pop; % Current Population
       PopSize = 2000; % size of the Population
       GenerationGenerator = 'NP'; % Method to generate next generation NORM (Parametric Normal Distribution), NP (Non Parametric Distribution), CR (Crossover and Mutation, traditional genetic algorithm)
       FacePenalty = 0;% Regularization term , the higher the number the more average the facial prediction (data versus model plausibility)
       ResFaceType = 'WeightedAverage';% WeightedAverage Average Best
       ResFace = [];% resulting face
       ResFaceC = [];% resulting face normalized coefficients
       ResFaceSpace = [];% If multiple evolution are performed, the output is a space stored in this property
       Scores = 0;% scores for all Pop individuals
       WeightGT = 'FREQ'; % Weighting the individual scores of matches to SNP models, NO, Atyp, FREQ
       InitFcn = 'Normal';% uniform in ranges, Gaussian according to eigenvalues,Training given populatio
       MaxGenerations = 100;
       Display = false;
       Verbose = false;
       EvolveCycles = 'Single';% the number of evolutions to perform: single or multiple
       nrCycles = 100;% Number of cycles to perform in case Evolve Cycles is set to multiple
       MatchUse = 0;% All (0) BF (1) SNPs (2)
    end
    properties (Dependent = true)
       Diversity;
       pFaces;
       BestScore;
       MeanScore;
       WorstScore;
       PopAverage; 
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
       SNPs;
       nrSNP;
       Dim;
       BestInd;
       WorstInd;
       sortIndex;
       GlobalRange;
       LocalRange;
       MatchW;
       GTW;
    end
    % REPRODUCTION
    properties 
       Reduction = 0.5; % Fraction to reduce
       ScalingFcn = 'Rank';% Scaling options
    end
    properties (Dependent = true)
       RedInd;
    end
    properties (Hidden = true)
       Tol = 10;
       mu = 0;
       Sigma = 0;
    end
    % LOCAL SOLVER
    properties 
       LocalSolver = 'fminsearch';
       SolveLocal = false;
       LocalSolution = [];
       LocalScore = 0;
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
        function obj = EvoMorphv2(varargin)
            obj = obj@superClass(varargin{:});
        end
    end
    methods % GENERAL GETTING/SETTING
        function out = get.SNPCont(obj)
           out = obj.SNPCont;
           if ~superClass.isH(out), out = []; end
        end
        function out = get.BaseCont(obj)
           out = obj.BaseCont;
           if ~superClass.isH(out), out = []; end
        end
        function out = get.nrSNP(obj)
            if isempty(obj.SNPCont), out = 0; return; end
            out = obj.SNPCont.nrSNP;
        end
        function out = get.PopSize(obj)
            out = obj.PopSize;
            if ~isempty(out), return; end
            out = 10*obj.Dim;% By default this is a good Pop size
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
        function out = get.SNPs(obj)
            if isempty(obj.SNPCont), out = []; return; end
            out = obj.SNPCont.SNPs;
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
        function out = get.GTW(obj)
            out = ones(obj.nrBF+obj.nrSNP,1);
            switch obj.WeightGT
                case 'ATYP'
                    out(obj.nrBF+1:end) = obj.GTAtyp;
                case 'FREQ'
                    out(obj.nrBF+1:end) = 1-obj.GTFreq;
                otherwise
                    return;
            end
            out(isnan(out)) = 0;
        end
        function out = get.MatchW(obj)
            out = obj.GTW.*obj.RW;
            if obj.MatchUse==0, return; end
            if obj.MatchUse==1, out = out(1:obj.nrBF); return; end
            out = out(obj.nrBF+1:end);
        end
        function out = get.nrCov(obj)
                 out = length(obj.Cov); 
        end
        function out = get.nrGB(obj)
                 out = length(obj.GB);
        end
        function out = get.nrBF(obj)
                 out = obj.nrCov+obj.nrGB;
        end
        function out = get.RedInd(obj)
                 out = ceil(obj.Reduction*obj.PopSize); 
        end
    end
    methods % MAIN INTERFACE FUNCTIONS
        function evolve(obj)
            switch obj.EvolveCycles
                case 'Single'
                    evolveSingle(obj);
                case 'Multiple'
                    obj.Display = false;
                    CyclePop = nan*zeros(obj.Dim,obj.nrCycles);
                    warning off;
                    tic;
                    %parfor_progress(obj.nrCycles);
                    FM = obj.FaceModel;
                    DM = obj.DemiModelSet;
                    parfor i=1:obj.nrCycles
                        out = [];
                        forobj = clone(obj);
                        forobj.FaceModel = FM;% passing the pointer not the whole object to save memory (see transient aspect)
                        forobj.DemiModelSet = DM;
                        forobj.EvolveCycles = 'Single';
                        forobj.Display = false;
                        evolve(forobj);
                        switch obj.ResFaceType
                            case 'WeightedAverage'
                                out = sum(repmat(forobj.scaledScores,1,forobj.Dim).*forobj.Pop)/sum(forobj.scaledScores);
                            case 'Average'
                                out = forobj.PopAverage;
                            case 'Best'
                                out = forobj.Pop(forobj.BestInd,:);
                        end
                        CyclePop(:,i) = out(:);
                        delete(forobj);
                        %parfor_progress;
                    end
                    toc;
                    warning on;
                    obj.PopSize = obj.nrCycles;
                    obj.Pop = CyclePop';
                    evalPop(obj);
                    scaleScores(obj);
                    getSolutionSpace(obj);
                    %parfor_progress(0);
            end
            updateResFace(obj);
        end
        function evolveSingle(obj)
            obj.State = 'Init';
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
            if obj.MovieRecord&&obj.Display
               close(obj.MovieFile);
            end
            if obj.SolveLocal,solveLocally(obj);end
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
            if strcmp(obj.State,'Init')
               obj.Scores = evalSolutions(obj,obj.Pop./repmat(obj.StdO',obj.PopSize,1),obj.pFaces);
            else
               red = obj.RedInd; 
               obj.Scores = obj.Scores(obj.sortIndex);
               redPop = obj.Pop(red+1:end,:);
               redpFaces = obj.pFaces(red+1:end);
               obj.Scores(red+1:end) = evalSolutions(obj,redPop./repmat(obj.StdO',size(redPop,1),1),redpFaces);
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
        function linkWithData(obj,basecont,snpcont,subject)
              if ~strcmp(basecont.TrackID,snpcont.TrackID), error('Base and SNP models are not compatible'); end
            % Is the required input present?
              obj.BaseCont = basecont;obj.SNPCont = snpcont;
            % Use Reduced Shape Space or NOT
              obj.UseRedShape = obj.SNPCont.SNPs{1}.ST4RedShape;
              if obj.UseRedShape
                   obj.FaceSpace = obj.BaseCont.RedShapeSpace;
                   obj.Cov = [];obj.GB = [];
              else
                   obj.FaceSpace = obj.BaseCont.ShapeSpace;
                   obj.Cov = subject.Cov;obj.GB = subject.BFGB;
                   [obj.CovRIP,obj.CovDistr] = CovX2RIP(obj.BaseCont,obj.Cov,subject.BFKappa,subject.BFK);
                   [obj.GBRIP,obj.GBDistr] = GBX2RIP(obj.BaseCont,obj.GB,subject.BFKappa,subject.BFK);
              end
             % Copy genotypes and other input parameters 
              obj.GT = subject.LinkedGT;
              obj.GTFreq = subject.GTFreq;
              obj.GTAtyp = subject.GTAtyp;
              obj.BF = subject.EvoBF;
              obj.ResFace = clone(obj.RefFace);
        end
        function out = evalSolutions(obj,sol,pSol)
            n = size(sol,1);
            if obj.MatchUse==0
                BaseMatches = -log(matchFaces(obj.BaseCont,sol,obj.CovDistr,obj.GBDistr));
                SNPMatches = -log(matchFaces(obj.SNPCont,sol,obj.GT,obj.BF));
                Matches = [BaseMatches;SNPMatches];
            elseif obj.MatchUse==1
                Matches = -log(matchFaces(obj.BaseCont,sol,obj.CovDistr,obj.GBDistr));
            else
                Matches = -log(matchFaces(obj.SNPCont,sol,obj.GT,obj.BF));
            end
            out = (1-obj.FacePenalty).*nansum(repmat(obj.MatchW,1,n).*Matches,1)./repmat(nansum(obj.MatchW),1,n)+obj.FacePenalty.*(pSol);
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
            W(obj.sortIndex(300:end)) = 0;
            out = sum(repmat(W,1,obj.Dim).*obj.Pop)/sum(W);
            %out = sum(repmat(obj.scaledScores,1,obj.Dim).*obj.Pop)/sum(obj.scaledScores);
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