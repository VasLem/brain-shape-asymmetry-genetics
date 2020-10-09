classdef JanusEvoMorph < superClass
    properties (Transient = true)% To make sure that certain object are not co-deleted, co-cloned or co-saved
        FaceModel;% Facial shape model, is not co-cloned nor deleted
        DemiModelSet;% List of Demi models
    end
    properties % overall properties
        Population; % Current Population
        PopulationSize = [];% Desired population size
        FacePenalty = 1;% Regularization term , the higher the number the more average the facial prediction (data versus model plausibility)
        ResFaceType = 'WeightedAverage';% WeightedAverage Average Best
        ResFace;% resulting face
        ResFaceSpace;% If multiple evolution are performed, the output is a space stored in this property
        GlobalScale = 2.5;% Global search space factor times standard deviations in face model
        LocalScale = 1;% Local search space factor times standard deviations in face model
        Scores;% scores for all population individuals
        InitFcn = 'Training';% uniform in ranges, Gaussian according to eigenvalues,Training given populatio
    end
    properties %re-production
        ElitePerc = 0.01; % Percentage of elite solutions to keep
        CrossoverFraction = 0.8;% Percentage of population-elite to produce using crossover, the rest is done using mutation
        CrossoverFcn = 'Scattered';% Scattered SinglePoint TwoPoint Intermediate Heuristic Arithmetic
        CrossoverIntermediateRatio = 0.5; % parameter for the Intermediate crossover function
        CrossoverHeuristicRatio = 1.2;% paramter for the heuristic crossover function
        MutationFcn = 'Normal';% Type of mutation function, options are Normal or Uniform
        MutationNormalShrink = 1; % Parameter to shrink the mutation space over generations (makes the population shrink)
        MutationUniformRate = 0.01;% Parameter when mutation is set to uniform
        SelectionFcn = 'StochUniform';% Selection of Parents: Uniform Remainder StochUniform Roulette Tournament
        TournamentSize = 4;% Tournament Parent selection
        ScalingFcn = 'Rank';% Rank Top Shift Proportional
        TopFrac = 0.4;% Top fraction to keep when using Top Scaling
        MaximumSurvivalRate = 2;
        BaseFaceFcn = 'X'; % Fcn to predict the base face, X = use the original predictor variables, RIP = take the step via RIP values
    end
    properties % algorithmic properties
        MaxGenerations = 100;
        MaxTime = +inf;
        FLimit = -inf;
        StallGenLimit = 50;
        %StallTimeLimit = Inf;
        FunctionTolerance = 0.0001;
        Display = false;
        vAxes = [1 2 3];
        TrimPopulation = false;
        LocalSolver = 'fminsearch';% fminsearch patternsearch
        SolveLocal = false;
        LocalSolution;% Currently this is to expensive, need to approximate objective function FIRST
        LocalScore;
        EvolveCycles = 'Single';% the number of evolutions to perform: single or multiple
        nrCycles = 100;% Number of cycles to perform in case Evolve Cycles is set to multiple
    end
    properties (Dependent = true)
       DMTags;
       DMLevels;
       nrDM;
       Diversity;
       pFaces;
       BestScore;
       MeanScore;
       WorstScore;
       PopulationAverage;
       Dim;
       DM;
       ID;
       Input;
       GlobalRange;
       LocalRange;
       nrParents;
       nrElite;
       nrCrossover;
       nrMutate;
    end
    properties (Dependent = true)% Prediction faces
        PredFace;
    end
    properties % Prediction faces and individuality assessments
        BaseFace;% Old school base face
        AmplFace;% amplified predicted face
        AF = 2;% amplification factor
        IdAss;% Identity assessment
    end
    properties (Hidden = true) % Algorithm parameters not to be seen by user
       Parents;
       scaledScores;
       OldBestScore = +inf;
       State = 'Init';
       Generation = 0;
       StallGeneration = 0;
       Time = 0;
       StartTime = 0;
       Fig;
       sError;
       sDiversity;
       sScores;
       sScaling;
       vFace;
       vPop;
       PopLM = [];
    end
    properties (Hidden = true, Dependent = true)
       BestInd;
       WorstInd;
       sortIndex;
       CrossoverParents;
       MutationParents;
       MutationScale;
       MutationRange;
       FaceSigmas;
    end
    methods % Constructor
        function obj = JanusEvoMorph(varargin)
            obj = obj@superClass(varargin{:});         
        end
    end
    methods % Special Setting & Getting
        function out = get.DMTags(obj)
           out = obj.DemiModelSet.DMTags; 
        end
        function out = get.DMLevels(obj)
           out = obj.DemiModelSet.DMLevels; 
        end
        function out = get.nrDM(obj)
            out = obj.DemiModelSet.nrUse;
        end
        function out = get.PopulationSize(obj)
            out = obj.PopulationSize;
            if ~isempty(out), return; end
            out = 10*obj.Dim;% By default this is a good population size
        end
        function out = get.Diversity(obj)
            out = mean(sqrt(sum((obj.Population-repmat(obj.PopulationAverage,obj.PopulationSize,1)).^2))); 
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
        function out = get.PopulationAverage(obj)
            out = mean(obj.Population);
        end
        function out = get.Dim(obj)
            if isempty(obj.FaceModel), out = 0; return; end
            out = obj.FaceModel.nrEV;
        end
        function out = get.FaceModel(obj)
           out = obj.FaceModel;
           if ~superClass.isH(out), out = []; end 
        end
        function out = get.DemiModelSet(obj)
           out = obj.DemiModelSet;
           if ~superClass.isH(out), out = []; end 
        end
        function out = get.ID(obj)
            if isempty(obj.DemiModelSet), out = []; return; end
            out = obj.DemiModelSet.ID;
        end
        function out = get.DM(obj)
            if isempty(obj.DemiModelSet), out = []; return; end
            out = obj.DemiModelSet.DM;
        end
        function obj = set.ID(obj,in)
            obj.DemiModelSet.ID = in;
        end
        function out = get.Input(obj)
            if isempty(obj.DemiModelSet), out = []; return; end
            out = obj.DemiModelSet.Input;
        end
        function obj = set.Input(obj,in)
            obj.DemiModelSet.Input = in;
        end
        function obj = set.DM(obj,in)
            obj.DemiModelSet.DM = in;
        end
        function out = get.GlobalRange(obj)
           if isempty(obj.FaceModel), out = []; return; end
           out = [obj.FaceModel.AvgCoeff-obj.GlobalScale*obj.FaceModel.EigStd obj.FaceModel.AvgCoeff+obj.GlobalScale*obj.FaceModel.EigStd]; 
        end
        function out = get.MutationRange(obj)
           if isempty(obj.FaceModel), out = []; return; end
           out = [obj.FaceModel.AvgCoeff-obj.MutationScale*obj.FaceModel.EigStd obj.FaceModel.AvgCoeff+obj.MutationScale*obj.FaceModel.EigStd]; 
        end
        function out = get.LocalRange(obj)
           if isempty(obj.FaceModel), out = []; return; end
           out = [obj.FaceModel.AvgCoeff-obj.LocalScale*obj.FaceModel.EigStd obj.FaceModel.AvgCoeff+obj.LocalScale*obj.FaceModel.EigStd]; 
        end
        function out = get.nrElite(obj)
           out = floor(obj.PopulationSize*obj.ElitePerc); 
        end
        function out = get.nrCrossover(obj)
            out = floor(obj.CrossoverFraction*(obj.PopulationSize-obj.nrElite));
        end
        function out = get.nrMutate(obj)
            out = obj.PopulationSize-obj.nrElite-obj.nrCrossover;
        end
        function out = get.nrParents(obj)
            out = obj.nrCrossover*2+obj.nrMutate;
        end
        function out = get.sortIndex(obj)
            if isempty(obj.Scores), out = []; end
            [~,out] = sort(obj.Scores);
        end
        function obj = set.FaceModel(obj,in)
            obj.FaceModel = in;
            obj.ResFace = clone(in.Average);
        end
        function out = get.CrossoverParents(obj)
           out = obj.Parents(1:obj.nrCrossover*2); 
        end
        function out = get.MutationParents(obj)
            out = obj.Parents(obj.nrCrossover*2+1:end);
        end
        function out = get.MutationScale(obj)
            out = obj.LocalScale - obj.MutationNormalShrink * obj.LocalScale * obj.Generation/obj.MaxGenerations;
        end
        function out = get.FaceSigmas(obj)
            out = obj.FaceModel.EigStd;
        end
        function out = get.PopLM(obj)
           out = obj.PopLM;
           if ~superClass.isH(out), out = []; end  
        end
        function out = get.pFaces(obj)
           out = sqrt(sum(((obj.Population-repmat(obj.FaceModel.AvgCoeff',obj.PopulationSize,1))./repmat(obj.FaceModel.EigStd',obj.PopulationSize,1)).^2,2))'; 
        end
        function out = get.PredFace(obj)
           updateResFace(obj); 
           out = obj.ResFace; 
        end
    end
    methods % General Interface Functions for the evolution proccess
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
                               out = sum(repmat(forobj.scaledScores,1,forobj.Dim).*forobj.Population)/sum(forobj.scaledScores);
                           case 'Average'
                               out = forobj.PopulationAverage;
                           case 'Best'
                               out = forobj.Population(forobj.BestInd,:);
                        end
                        CyclePop(:,i) = out(:);
                        delete(forobj);
                        %parfor_progress;
                    end
                    toc;
                    warning on;
                    obj.PopulationSize = obj.nrCycles;
                    obj.Population = CyclePop';
                    evalPopulation(obj);
                    scaleScores(obj);
                    getSolutionSpace(obj);
                    %parfor_progress(0);
            end
            updateResFace(obj);
        end
        function evolveSingle(obj)
            obj.State = 'Init';
            while ~strcmp(obj.State,'Stop');
                if strcmp(obj.State,'Init')
                   obj.StartTime = cputime;
                   obj.Time = 0;
                   obj.Generation = 1; 
                   InitializePopulation(obj);            
                else
                   obj.Generation = obj.Generation + 1;
                   obj.Time = cputime-obj.StartTime;
                   nextGeneration(obj);
                end
                evalPopulation(obj);
                scaleScores(obj);
                selectParents(obj);
                if obj.Display, 
                   plotEvolution(obj);
                else
                   %disp(['G: ' num2str(obj.Generation) ' BS: ' num2str(obj.BestScore) ' MS: ' num2str(obj.MeanScore)]);
                end
                stopTest(obj);
            end
            if obj.SolveLocal
               solveLocally(obj); 
            end
        end
        function solveLocally(obj)
            switch obj.ResFaceType
                case 'WeightedAverage'
                    X0 = sum(repmat(obj.scaledScores,1,obj.Dim).*obj.Population)/sum(obj.scaledScores);
                case 'Average'
                    X0 = obj.PopulationAverage;
                case 'Best'
                    X0 = obj.Population(obj.BestInd,:);
            end
            options = optimset('Display','Iter');
            switch obj.LocalSolver
                case 'fminsearch'
                    [obj.LocalSolution,obj.LocalScore] = fminsearch(@(X) myLocalSolverFun(obj,X),X0,options);
                case 'patternsearch'
                    [obj.LocalSolution,obj.LocalScore] = patternsearch(@(X) myLocalSolverFun(obj,X),X0,[],[],[],[],[],[],[],options);
            end
        end
        function InitializePopulation(obj)
            if isempty(obj.FaceModel), error('First Load Face Model'); end
            switch obj.InitFcn
                case 'Uniform'
                    R = obj.GlobalRange;
                    obj.Population = (repmat(R(:,1),1,obj.PopulationSize) + repmat(R(:,2)-R(:,1),1,obj.PopulationSize).*rand(obj.Dim,obj.PopulationSize))';
                case 'Normal'
                    obj.Population = (repmat(obj.FaceModel.AvgCoeff,1,obj.PopulationSize) + repmat(obj.FaceModel.EigStd,1,obj.PopulationSize).*randn(obj.Dim,obj.PopulationSize))';
                case 'Training'
                    ind = randsample(obj.FaceModel.n,obj.PopulationSize,true);
                    obj.Population = obj.FaceModel.Tcoeff(ind,:);
            end
        end    
        function evalPopulation(obj)
            obj.Scores = getScore(obj.DemiModelSet,obj.Population') + obj.FacePenalty*obj.pFaces;
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
        function selectParents(obj)
            switch obj.SelectionFcn
                case 'Uniform'
                    obj.Parents = selectionuniform(obj); 
                case 'Remainder'
                    obj.Parents = selectionremainder(obj);
                case 'StochUniform'
                    obj.Parents = selectionstochunif(obj);
                case 'Roulette'
                    obj.Parents = selectionroulette(obj);
                case 'Tournament'
                    obj.Parents = selectiontournament(obj);
            end  
        end
        function nextGeneration(obj)
            % elite selection
            if obj.nrElite>0
               elite = obj.Population(obj.selectElite,:);
            else
               elite = [];
            end
            % crossover children
            if obj.nrCrossover>0
               Xover = crossoverChildren(obj); 
            else
               Xover = [];
            end
            % mutation children
            if obj.nrMutate>0
               Mut = mutateChildren(obj); 
            else
               Mut = [];
            end
            obj.Population = [elite; Xover; Mut];
            trimPopulation(obj);
        end
        function out = selectElite(obj)
            out = obj.sortIndex(1:obj.nrElite);
        end
        function out = crossoverChildren(obj)
            switch obj.CrossoverFcn
                case 'Scattered'
                    out  = crossoverscattered(obj);
                case 'SinglePoint'
                    out = crossoversinglepoint(obj);
                case 'TwoPoint'
                    out = crossovertwopoint(obj);
                case 'Intermediate'
                    out = crossoverintermediate(obj);
                case 'Heuristic'
                    out = crossoverheuristic(obj);
                case 'Arithmetic'
                    out = crossoverarithmetic(obj);
            end
        end
        function out = mutateChildren(obj)
            switch obj.MutationFcn
                case 'Normal'
                    out = mutationgaussian(obj);
                case 'Uniform'
                    out = mutationuniform(obj);
            end
        end
        function stopTest(obj)
            obj.State = 'Run';
            if obj.Generation>obj.MaxGenerations, obj.State = 'Stop'; end
            if obj.Time>obj.MaxTime, obj.State = 'Stop'; end
        end
        function trimPopulation(obj)
            % population must remain within global range
            if ~obj.TrimPopulation, return; end
            M = repmat(obj.GlobalRange(:,1)',obj.PopulationSize,1);
            index = find(obj.Population<M);
            if ~isempty(index),obj.Population(index) = M(index);end
            M = repmat(obj.GlobalRange(:,2)',obj.PopulationSize,1);
            index = find(obj.Population>M);
            if ~isempty(index),obj.Population(index) = M(index);end
        end
        function out = feasibleKid(obj,in)
            out = true;% To be implemented in case this is needed
        end
    end
    methods % Plotting functions 
        function updateResFace(obj)
           switch obj.ResFaceType
               case 'WeightedAverage'
                   out = getWAverageFace(obj);
               case 'Average'
                   out = getAverageFace(obj);
               case 'Best'
                   out = getBestFace(obj);
           end
           obj.ResFace.Vertices = out.Vertices;
           if strcmp(obj.FaceModel.Type,'appearancePCA')
               obj.ResFace.TextureMap = clone(out.TextureMap);
           end
        end   
        function updatePopLM(obj)
            obj.PopLM.Vertices = obj.Population(:,obj.vAxes)';
            obj.PopLM.Value = obj.Scores;
        end
        function plotEvolution(obj)
            if strcmp(obj.State,'Init');
               constructDisplay(obj);
            end
            updateResFace(obj);
            updatePopLM(obj);
            plot(obj.sError,obj.Generation,obj.BestScore,'k.');
            plot(obj.sError,obj.Generation,obj.MeanScore,'b.');
            title(obj.sError,['Best score: ' num2str(obj.BestScore) ' Mean Score: ' num2str(obj.MeanScore)]);
            plot(obj.sDiversity,obj.Generation,obj.Diversity,'b.');
            hist(obj.sScores,obj.Scores,20);
            plot(obj.sScaling,obj.Scores,obj.scaledScores,'b.');
            drawnow;
        end
    end
    methods % Extracting range of facial solutions   
        function out = getWorstFace(obj)
            out = getScan(obj.FaceModel,obj.Population(obj.WorstInd,:));
        end
        function out = getAverageFace(obj)
            out = getScan(obj.FaceModel,obj.PopulationAverage);
        end
        function out = getBestFace(obj)
            out = getScan(obj.FaceModel,obj.Population(obj.BestInd,:));
        end
        function out = getLocalFace(obj)
            out = getScan(obj.FaceModel,obj.LocalSolution);
        end
        function out = getWAverageFace(obj)
            out = getWAverage(obj);
            out = getScan(obj.FaceModel,out);
        end
        function out = getWAverage(obj)
            out = sum(repmat(obj.scaledScores,1,obj.Dim).*obj.Population)/sum(obj.scaledScores);
        end
        function out = getSolutionSpace(obj)
           out = clone(obj.FaceModel);
           out.Tcoeff = obj.Population;
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
        function out = getBaseFace(obj,B)
           % this function is a reimplementation of the two-stage prediction strategy as published in FSI-genetics
           % The base face is the face without the effect of individual
           % SNPs, I have opted to implement this ability within the
           % JANUSEvoMorph class indespite it has notthing to do with
           % Evolutionary computing, however it is a valuable tool in the
           % assesment of the individuality of the prediction
           if nargin < 2, B = obj.FaceModel.Tcoeff; end
           nrTr = size(B,1);
           DMindex = find(strcmp('Genomic',obj.DMLevels));
           nrDM = length(DMindex);
           TrA = zeros(nrTr,nrDM);
           A = zeros(1,nrDM);
           switch obj.BaseFaceFcn
               case 'X' % simple one step regression
                   for i=1:1:nrDM
                       TrA(:,i) = obj.DemiModelSet.DM{DMindex(i)}.TrX;
                       A(i) = obj.Input(DMindex(i));
                   end
               case 'RIP' % more complicated regression via RIP values
                   for i=1:1:nrDM
                       DM = obj.DemiModelSet.DM{DMindex(i)};
                       TrA(:,i) = DM.EstRIP;
                       switch DM.XType
                           case 'Categorical'
                               A(i) = DM.CatMu(DM.CatLabel==obj.Input(DMindex(i)));
                           case 'Continous'
                               [X,RIP] = eliminateNAN(DM.TrX,DM.EstRIP');
                               P = polyfit(X,RIP,DM.D);
                               A(i) = polyval(P,obj.Input(DMindex(i)));
                       end
                   end
           end
           % from now on, simply perform a multiple regression and apply it
           % to the desired values
           [TrA,B] = eliminateNAN(TrA,B);
           [~,~,~,~,MC] = plsregress(TrA,B,size(TrA,2));
           mC = mean(TrA);
           mShape = mean(B);
           DX = mC-A;
           DY = DX*MC(2:end,:);
           NY = mShape-DY;
           out = getScan(obj.FaceModel,NY');
           obj.BaseFace = out;
        end
    end
    methods % Other general interface functions
        function parseInput(obj,in)
           if isempty(obj.DemiModelSet), error('First Load DemiModelSet'); end
           parseInput(obj.DemiModelSet,in); 
        end
        function delete(obj)
            try 
                close(obj.Fig);
            catch
            end
            %superClass.delete(obj);
        end
    end
    methods % Main Prediction function + amplifying and exporting functions
       function getFacialPrediction(obj)
                getBaseFace(obj);
                switch obj.EvolveCycles
                    case 'Single'
                        evolve(obj);
                    otherwise 'Multiple'
                        out = 1;
                end
       end 
       function amplifyIdentity(obj)
                disp('Amplifying Identity...');
                if isempty(obj.PredFace), error('Make a prediction first'); end
                dir = obj.PredFace.Vertices-obj.BaseFace.Vertices;
                obj.AmplFace = clone(obj.PredFace);
                obj.AmplFace.Vertices = obj.BaseFace.Vertices + obj.AF*dir;
                disp('--------------------------------------------------------------');
       end
       function assessIdentity(obj)
                disp('Assessing Identity...');
                obj.IdAss = compareMorphs(obj.AmplFace,obj.BaseFace);
                disp('--------------------------------------------------------------');
       end
       function exportFaces(obj)
                disp('Exporting Faces..');
                exportWavefront(obj.BaseFace,'BaseFace.obj');
                exportWavefront(obj.PredFace,'PredictedFace.obj');
                exportWavefront(obj.AmplFace,'AmplifiedFace.obj');
                disp('--------------------------------------------------------------');
       end
    end
end