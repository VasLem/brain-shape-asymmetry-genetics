classdef G2FPredictionValidation < superClass
   properties
      ID = [];
      S = nan;
      A = nan;
      W = nan;
      H = nan;
      GB = [];
      GT = [];
      GBNames = {};
      RS = {};
      BaseModel = 'EURO';
      SNPModel = 'EURO';
      BaseFace = [];
      PredFace = [];
      AmplFace = [];
      SNPUsed = [];
   end
   properties (Dependent = true)
      nSNPUsed; 
      AvgOcc;
      AvgAtyp;
      nrMAtyp;
   end
   properties (Hidden = true)
       nGB = [];
       nSNP = [];
   end
   properties (Dependent = true, Hidden = true)
      Cov;
      SRIP;
      ARIP;
      HRIP;
      WRIP;
   end
   % BASE FACE CREATION
   properties
      BFType = 'RIP'; % RIP or X
      BFAdjust = 'RESCALE';% RESCALE/ TRIM, NO
   end
   properties (Hidden = true)
      BFKappa = 3;
      BFK = 300;
      BFC = [];
      BFAtyp = [];
      BFTrim = 2;
   end
   properties (Transient = true, Hidden = true) 
      BFNrGB = 0;
      BFGBExtract = [];
      BFGBPresent = [];
      BFCovRIP = [];
      BFCovDistr = {};
      BFGBRIP = [];
      BFGBDistr = {};
      BFTrackID = [];
   end
   properties (Dependent = true, Hidden = true)
      BFGB; 
   end
   % GENOTYPE PREDICTIONS
   properties 
      PredMethod = 'OVERLAYING';% OVERLAYING/ EVO 
   end
   properties (Transient = true, Hidden = true)
      NrSNP = [];
      SNPExtract = [];
      SNPPresent = [];
      GTRIP = [];
      GTDistr = [];
      SNPContTrackID = [];
      SNPRedShape = false;
      SNPRef = [];
   end
   properties (Hidden = true)
      GTInfo = [];
      GTAtyp = []; 
      SNPPredTrim = 2;
      PredC = [];
   end
   properties (Dependent = true, Hidden = true)
      LinkedGT; 
   end
   % SNP OVERLAYING
   properties
      SOType = 'SEQ';% AVG (Averaging) or SEQ (Sequential) 
      SOAdjust = 'RESCALE';% NO/ RESCALE or TRIM
   end
   properties (Hidden = true)
      SOBF = 0;
      SOMAtyp = 2;
      SOSEQRuns = 10;
      SOW = 'NO';%LM (R2 per landmark)/CO (R2 per landmark coordinate)/ATYP (morph atyp value)/PC (R2 per PC)/NO (no weighting)  
   end
   % SNP EVOLUTION
   properties (Transient = true, Hidden = true)
      Evo = [];
   end
   properties (Hidden = true)
      EvoBF = 0;
   end
   % TEXTURE EXTRACTION
   properties
       TexMethod = 'BF';% BF using BaseFace Information, BFSHAPE , texture from Base Shape, PREDSHAPE, texture from Predicted Face
   end
   % IDENTITY ASSESSMENT AND AMPLIFICATION
   properties
      AF = 2;
      IdAss = [];
   end
   methods % Constructor
        function obj = G2FPredictionValidation(varargin)
            obj = obj@superClass(varargin{:});
            obj.Evo = EvoMorph;
        end
   end
   methods % GENERAL GETTING\SETTING
        function out = get.nSNPUsed(obj)
           out = length(obj.SNPUsed);
        end
        function out = get.Evo(obj)
            out = obj.Evo;
            if ~superClass.isH(out), out = []; end
        end
        function out = get.AvgOcc(obj)
            if isempty(obj.GTInfo), out = nan; return; end
            out = nanmean(obj.GTInfo.Occ);
        end
        function out = get.AvgAtyp(obj)
            if isempty(obj.GTInfo), out = nan; return; end
            out = nanmean(obj.GTInfo.atyp);
        end
        function out = get.nrMAtyp(obj)
            if isempty(obj.GTInfo), out = nan; return; end
            out = (nansum(obj.GTInfo.matyp)/obj.nSNPUsed)*100;
        end
        function out = get.Cov(obj)
                 out = [obj.S obj.A obj.W obj.H]; 
        end
        function out = get.SRIP(obj)
            if isempty(obj.CovRIP), out = []; return; end
            out = obj.CovRIP(1);
        end
        function out = get.ARIP(obj)
            if isempty(obj.CovRIP), out = []; return; end
            out = obj.CovRIP(2);
        end
        function out = get.WRIP(obj)
            if isempty(obj.CovRIP), out = []; return; end
            out = obj.CovRIP(3);
        end
        function out = get.HRIP(obj)
            if isempty(obj.CovRIP), out = []; return; end
            out = obj.CovRIP(4);
        end
   end
   methods % BASE FACE GETTING\SETTING
       function out = get.BFGB(obj)
          out = [];  
          if isempty(obj.BFNrGB), return; end
          out = nan*zeros(1,obj.BFNrGB);
          if isempty(obj.BFGBPresent), return; end
          if isempty(obj.BFGBExtract), return; end
          out(obj.BFGBPresent) = obj.GB(obj.BFGBExtract);
       end
   end
   methods % SNP OVERLAYING GETTING\SETTING
       function out = get.LinkedGT(obj)
          out = []; 
          if isempty(obj.NrSNP), return; end
          if isempty(obj.SNPPresent), return; end
          if isempty(obj.SNPExtract), return; end
          out = nan*zeros(1,obj.NrSNP);
          out(obj.SNPPresent) = obj.GT(obj.SNPExtract);
       end
   end
   methods % MAIN INTERFACE FUNCTIONS
       function FaceYourDNA(obj,BaseCont,SNPCont,TexCont,SNPBaseCont)
                if nargin < 5, SNPBaseCont = BaseCont; end
                % STEP 1: BASE FACE
                  linkWithBaseCont(obj,BaseCont);
                  extractBaseFace(obj,BaseCont);
                % STEP 2: PREDICTED FACE
                  linkWithSNPCont(obj,SNPCont);
                  extractPredictedFace(obj,SNPBaseCont,SNPCont);
                % STEP 3: TEXTURE
                  extractTexture(obj,TexCont,BaseCont,SNPCont);
       end
       function parsePredictionInput(obj,input)
           disp('Parsing Prediction Input...');
           if datenum(date)>datenum('31-Oct-2015'), error('License Expired, please contact peter.claes@esat.kuleuven.be');end
           if datenum(date)+30>datenum('31-Oct-2015')
              warning(['License will expire in ' num2str(datenum('31-Oct-2015')-datenum(date)) ' days , please contact peter.claes@esat.kuleuven.be']);
           end
           [~,~,raw] = xlsread(input);
           obj.nGB = raw{1,2};
           obj.nSNP = raw{2,2};
           %obj.WithRIPG = logical(raw{3,2});
           obj.ID = raw{4,2};
           obj.S = raw{5,2};
           obj.A = raw{6,2};
           obj.H = raw{7,2};
           obj.W = raw{8,2};
           tmp = raw(9:8+obj.nGB,2)';
           obj.GB = nan(size(tmp));
           for i=1:1:length(tmp)
               if isnumeric(tmp{i})
                  obj.GB(i) = tmp{i};
               end
           end    
           obj.GBNames = raw(9:8+obj.nGB,1)';
           obj.RS = raw(9+obj.nGB:8+obj.nGB+obj.nSNP,1)';
           tmp = raw(9+obj.nGB:8+obj.nGB+obj.nSNP,2)';
           obj.GT = nan(size(tmp));
           for i=1:1:length(tmp)
               if isnumeric(tmp{i})
                  obj.GT(i) = tmp{i};
               end
           end
           disp('--------------------------------------------------------------');
       end
       function extractBaseFace(obj,BaseCont)
                disp('Creating Base Face...');
                switch obj.BFType
                    case 'X'
                        [obj.BaseFace,obj.BFC,obj.BFAtyp] = getBaseFace(BaseCont,obj.Cov,obj.BFGB,'X',obj.BFAdjust,obj.BFTrim);
                    case 'RIP'
                        [obj.BaseFace,obj.BFC,obj.BFAtyp] = getBaseFace(BaseCont,obj.BFCovRIP,obj.BFGBRIP,'RIP',obj.BFAdjust,obj.BFTrim);
                end
                obj.BaseFace.SingleColor = [0.4 0.4 0.4];obj.BaseFace.Material = 'Dull';
                disp('----------------------------------------------------------------------------');
       end
       function extractPredictedFace(obj,BaseCont,SNPCont)
                switch lower(obj.PredMethod)
                    case 'overlaying'
                        overlaySNPs(obj,BaseCont,SNPCont);
                    case 'evo'
                        evolveSNPs(obj,BaseCont,SNPCont);
                end
                obj.PredFace.Value = vDistances(obj.PredFace,obj.BaseFace);
                obj.PredFace.SingleColor = [0.7 0.7 0.7];obj.PredFace.Material = 'Dull';
       end    
       function extractTexture(obj,TexCont,BaseCont,SNPCont)
                if ~strcmp(BaseCont.TrackID,TexCont.TrackID), error('Base Model and Texture Model are NOT compatible'); end
                switch lower(obj.TexMethod)
                    case 'bf'             
                        switch obj.BFType
                            case 'RIP'
                                TRX = [BaseCont.CovRIP, BaseCont.GBRIP];
                                X = [obj.BFCovRIP, obj.BFGBRIP];
                            otherwise
                                TRX = [BaseCont.Cov BaseCont.GB];
                                X = [obj.Cov obj.BFGB];
                        end
                    case 'bfshape'
                        TRX = [BaseCont.DepVar];
                        X = obj.BFC;
                    case 'predshape'
                        if ~strcmp(SNPCont.TrackID,TexCont.TrackID), error('SNP Model and Texture Model are NOT compatible'); end
                        TRX = [BaseCont.DepVar];
                        X = obj.PredC;   
                end
                out = extractTexture(TexCont,TRX,X);
                obj.PredFace.UV = TexCont.UV;
                obj.PredFace.TextureMap = clone(out.TextureMap);
                obj.PredFace.ColorMode = 'Texture';
       end
   end
   methods % SUB INTERFACE FUNCTIONS
       function linkWithBaseCont(obj,BaseCont)
                disp('LINKING BASE MODEL');
                obj.BFNrGB = BaseCont.nrGB;
                [obj.BFGBExtract,obj.BFGBPresent] = vlookup(BaseCont.GBNames,obj.GBNames,false);
                [obj.BFCovRIP,obj.BFCovDistr] = CovX2RIP(BaseCont,obj.Cov,obj.BFKappa,obj.BFK);
                [obj.BFGBRIP,obj.BFGBDistr] = GBX2RIP(BaseCont,obj.BFGB,obj.BFKappa,obj.BFK);
                obj.BFTrackID = BaseCont.TrackID;
                disp('----------------------------------------------------------------------------');
       end
       function linkWithSNPCont(obj,SNPCont)
                disp('LINKING SNP MODEL');
                obj.NrSNP = SNPCont.nrSNP;
                [obj.SNPExtract,obj.SNPPresent] = vlookup(getRS(SNPCont),obj.RS,false);
                [obj.GTRIP,obj.GTDistr] = GT2RIP(SNPCont,obj.LinkedGT);
                obj.SNPContTrackID = SNPCont.TrackID;
                obj.GTAtyp = getGTAtyp(SNPCont,obj.LinkedGT);
                disp('----------------------------------------------------------------------------');
       end
       function overlaySNPs(obj,BaseCont,SNPCont)
                if isempty(obj.BaseFace), error('First create BaseFace or call getFacialPrediction');end
                if ~strcmp(BaseCont.TrackID,SNPCont.TrackID), error('Base Model and SNP Model are NOT compatible'); end
                disp('Matching Genotypes & overlaying effects...');
                switch lower(obj.SOType)
                    case {'avg' 'add'}
                        out = overlayAvgAddSNPs(obj,BaseCont,SNPCont);
                    case 'seq'
                        out = overlaySeqSNPs(obj,BaseCont,SNPCont);
                end
                obj.PredFace = clone(obj.BaseFace);
                obj.PredFace.Vertices = obj.PredFace.Vertices + out.ShapeRes;
                obj.PredFace.Value = vDistances(obj.PredFace,obj.BaseFace);
                obj.SNPUsed = find(out.Good);
                obj.PredC = out.C;
                disp([num2str(obj.nSNPUsed) ' SNPs Used']);
                disp('--------------------------------------------------------------');
       end
       function out = overlaySeqSNPs(obj,BaseCont,SNPCont)
                good = zeros(1,SNPCont.nrSNP);
                obj.GTInfo.Occ = nan*zeros(1,SNPCont.nrSNP);
                obj.GTInfo.atyp = nan*zeros(1,SNPCont.nrSNP);
                obj.GTInfo.matyp = nan*zeros(1,SNPCont.nrSNP);
                obj.GTInfo.MorphAtyp = nan*zeros(1,SNPCont.nrSNP);
                genotypes = obj.LinkedGT;
                % get order of SNP strength
                PS = getValues(SNPCont,'ST3PStrength');
                [~,sortindex] = sort(PS,'ascend');
                % Gather starting information by selecting one of the snp models
                obj.SNPRedShape = SNPCont.SNPs{1}.ST4RedShape;
                switch obj.SNPRedShape
                       case 1
                          Ref.Coeff = BaseCont.ST4Ref.RedCoeff;
                          Ref.Scan = clone(BaseCont.ST4Ref.RedScan);
                          Ref.Cov = [];
                       case 0
                          Ref.Coeff = BaseCont.ST4Ref.Coeff;
                          Ref.Scan = clone(BaseCont.ST4Ref.Scan);
                          switch BaseCont.RedShapeSpaceType
                              case 'X'
                                  Ref.Cov = [BaseCont.ST4Ref.AvgCov, BaseCont.ST4Ref.AvgGB];
                              case 'RIP'
                                  Ref.Cov = [BaseCont.ST4Ref.CovRIP, BaseCont.ST4Ref.GBRIP];
                          end
                end
                InitRef = Ref;
                f = waitbar(0,'Processing...');drawnow;
                for i=1:1:obj.SOSEQRuns
                    for g=1:1:SNPCont.nrSNP
                        gt = genotypes(sortindex(g));if isnan(gt), continue; end
                        [C] = overLaySNP(SNPCont.SNPs{sortindex(g)},BaseCont,gt,Ref,obj.SOBF,obj.SOMAtyp,'NO',10);
                        if isempty(C), continue; end
                        [obj.GTInfo.Occ(sortindex(g)),obj.GTInfo.atyp(sortindex(g)),obj.GTInfo.matyp(sortindex(g))] = GTTypicality(SNPCont.SNPs{sortindex(g)},gt);
                        Ref.Coeff = C;
                        good(sortindex(g)) = 1;
                    end
                    waitbar(i/obj.SOSEQRuns,f);drawnow;
                end
                if obj.SNPRedShape
                   ShapeSpace = BaseCont.RedShapeSpace;
                else
                   ShapeSpace = BaseCont.ShapeSpace;
                end
                out.C = C;
                switch lower(obj.SOAdjust);
                    case 'rescale'
                        C = C.*ShapeSpace.EigStd';
                    case 'trim'
                        B = obj.SNPPredTrim*ShapeSpace.EigStd';
                        index = find(abs(C)>B);
                        C(index) = B(index).*sign(C(index));
                    otherwise
                end
                %out.C = C;
                out.Scan = getScan(ShapeSpace,C);
                out.ShapeRes = out.Scan.Vertices-InitRef.Scan.Vertices;
                out.Good = good;
                delete(f);
       end
       function out = overlayAvgAddSNPs(obj,BaseCont,SNPCont)
                f = waitbar(0,'Processing...');drawnow;
                good = zeros(1,SNPCont.nrSNP);
                shapeR = nan*zeros(obj.BaseFace.nrV*3,SNPCont.nrSNP);
                weightsR = nan*zeros(obj.BaseFace.nrV*3,SNPCont.nrSNP);
                obj.GTInfo.Occ = nan*zeros(1,SNPCont.nrSNP);
                obj.GTInfo.atyp = nan*zeros(1,SNPCont.nrSNP);
                obj.GTInfo.matyp = nan*zeros(1,SNPCont.nrSNP);
                obj.GTInfo.MorphAtyp = nan*zeros(1,SNPCont.nrSNP);
                genotypes = obj.LinkedGT;
                for g=1:1:SNPCont.nrSNP
                    gt = genotypes(g);if isnan(gt), continue; end
                    if g==1 % Figure out reference to use
                       obj.SNPRedShape = SNPCont.SNPs{g}.ST4RedShape;
                       switch obj.SNPRedShape
                           case 1
                               Ref.Coeff = BaseCont.ST4Ref.RedCoeff;
                               Ref.Scan = clone(BaseCont.ST4Ref.RedScan);
                               Ref.Cov = [];
                           case 0
                               Ref.Coeff = BaseCont.ST4Ref.Coeff;
                               Ref.Scan = clone(BaseCont.ST4Ref.Scan);
                               switch BaseCont.RedShapeSpaceType
                                    case 'X'
                                       Ref.Cov = [BaseCont.ST4Ref.AvgCov, BaseCont.ST4Ref.AvgGB];
                                    case 'RIP'
                                       Ref.Cov = [BaseCont.ST4Ref.CovRIP, BaseCont.ST4Ref.GBRIP];
                               end
                       end
                    end
                    [~,Atyp,morph] = overLaySNP(SNPCont.SNPs{g},BaseCont,gt,Ref,obj.SOBF,obj.SOMAtyp,obj.SOAdjust,obj.SNPPredTrim);
                    if isempty(morph), continue; end
                    [obj.GTInfo.Occ(g),obj.GTInfo.atyp(g),obj.GTInfo.matyp(g)] = GTTypicality(SNPCont.SNPs{g},gt);
                    %shapeR(:,g) = morph.Vertices(:)-Ref.Scan.Vertices(:)+obj.BaseFace.Vertices(:);
                    shapeR(:,g) = morph.Vertices(:)-Ref.Scan.Vertices(:);
                    switch obj.SOW
                        case 'LM'
                            tmp = repmat(SNPCont.SNPs{g}.ST4MOD.R2LM,3,1);
                            weightsR(:,g) = tmp(:);
                        case 'CO'
                            weightsR(:,g) = SNPCont.SNPs{g}.ST4MOD.R2(:);
                        case 'NO'
                            weightsR(:,g) = ones(size(weightsR(:,g)));
                        case 'ATYP'
                            weightsR(:,g) = Atyp*ones(size(weightsR(:,g)));
                        otherwise
                    end
                    good(g) = 1;
                    waitbar(g/SNPCont.nrSNP,f);drawnow;
                end
                shapeR = shapeR(:,good==1);
                weightsR = weightsR(:,good==1);
                switch lower(obj.SOType)
                    case 'avg'
                        disp('Averaging');
                        shaperes = sum(shapeR.*weightsR,2)./sum(weightsR,2);
                    case 'add'
                        disp('Additive');
                        shaperes = sum(shapeR.*weightsR,2);
                end
                out.ShapeRes = reshape(shaperes,3,Ref.Scan.nrV);
                out.Scan = clone(Ref.Scan);
                out.Scan.Vertices = out.Scan.Vertices+out.ShapeRes;
                out.Good = good;
                out.C = [];
                delete(f);
       end
       function evolveSNPs(obj,BaseCont,SNPCont)
                obj.Evo.Display = true;
                obj.Evo.PopSize = 4000;
                obj.Evo.GenerationGenerator = 'NP';% 'NP'
                obj.Evo.FacePenalty = 0;
                obj.Evo.WeightMatches = 'NO';
                obj.Evo.MaxGenerations = 150;
                obj.Evo.ResFaceType = 'WeightedAverage';
                %obj.Evo.ResFaceType = 'Best';
                obj.Evo.GlobalScale = 3;
                obj.Evo.Reduction = 0.6;
                obj.Evo.InitFcn = 'Normal';
                linkWithData(obj.Evo,BaseCont,SNPCont,obj);
                evolve(obj.Evo);
                obj.PredFace = clone(obj.BaseFace);
                scan = getBestFace(obj.Evo);
                obj.PredFace.Vertices = scan.Vertices;
                obj.PredFace.Vertices = obj.Evo.ResFace.Vertices;
                obj.PredC = obj.Evo.ResFaceC;
                obj.PredFace.Value = vDistances(obj.PredFace,obj.BaseFace);
                obj.PredFace.SingleColor = [0.7 0.7 0.7];obj.PredFace.Material = 'Dull';
       end
   end
   methods % IDENTITY ASSESSMENT
       function amplifyIdentity(obj)
                disp('Amplifying Identity...');
                if isempty(obj.PredFace), error('Make a prediction first'); end
                dir = obj.PredFace.Vertices-obj.BaseFace.Vertices;
                obj.AmplFace = clone(obj.PredFace);
                obj.AmplFace.Vertices = obj.BaseFace.Vertices + obj.AF*dir;
                obj.AmplFace.SingleColor = [0.9 0.9 0.9];
                obj.AmplFace.Value = vDistances(obj.AmplFace,obj.BaseFace);
                disp('--------------------------------------------------------------');
       end
       function assessIdentity(obj)
                disp('Assessing Identity...');
                obj.IdAss = compareMorphs(obj.AmplFace,obj.BaseFace);
                disp('--------------------------------------------------------------');
                
       end
   end
   methods % IMAGING AND EXPORTING
       function imageFacePrediction(obj,RefScanV,vref,name,saveresults)
            mapname = 'jet';
            BaseFace = obj.BaseFace; %#ok<*PROP>
            AmplFace = obj.AmplFace;
            comparison = obj.IdAss;
            BaseFace.Material = 'Dull';
            AmplFace.Material = 'Dull';
%             % imaging the scan
%             [v,RefScan] = G2FPredictionV2.cloneAndDisplay(vref,RefScanV,'Vertices',BaseFace.Vertices,'Texture',[name '_Base'],[]);
%             if saveresults,G2FPredictionV2.savingResults(v,RefScan);end
            % imaging the norm
            [v,RefScan] = G2FPredictionV2.cloneAndDisplay(vref,RefScanV,'Vertices',AmplFace.Vertices,'Texture',[name '_Pred'],[]);
            if saveresults,G2FPredictionV2.savingResults(v,RefScan);end
            % imaging the magnitude of the displacements
%             val = vDistances(BaseFace,AmplFace);
%             [v,RefScan] = G2FPredictionV2.cloneAndDisplay(vref,AmplFace,'Value',val,'Indexed',[name 'Disp (Blue NO disp; Red MAX disp)'],[0 max(val)]);
%             colormap(mapname);
%             if saveresults,G2FPredictionV2.savingResults(v,RefScan);end
            % imaging the displacements along the normals
            [v,RefScan] = G2FPredictionV2.cloneAndDisplay(vref,RefScanV,'Value',comparison.NormalDistances,'Indexed',[name ' Normal (Blue Inward disp; Red Outward disp)'],[]);
            colormap(mapname);
            if saveresults,G2FPredictionV2.savingResults(v,RefScan);end
            % imaging the Area Ratios
            [v,RefScan] = G2FPredictionV2.cloneAndDisplay(vref,RefScanV,'Value',comparison.AreaRatios,'Indexed',[name ' Area (Blue Decrease; Red Increase )'],[]);
            colormap(mapname);
            if saveresults,G2FPredictionV2.savingResults(v,RefScan);end
            % imaging the signed Curvatures
            maxcurv = min(max(abs(comparison.signedCurvature1)),max(abs(comparison.signedCurvature2)));
%             [v,RefScan] = G2FPredictionV2.cloneAndDisplay(vref,BaseFace,'Value',comparison.signedCurvature2,'Indexed',[name ' Sign. Curv. Base (Blue Concave; Red Convex)'],[-1*maxcurv maxcurv]);
%             colormap(mapname);
%             if saveresults,G2FPredictionV2.savingResults(v,RefScan);end
            [v,RefScan] = G2FPredictionV2.cloneAndDisplay(vref,AmplFace,'Value',comparison.signedCurvature1,'Indexed',[name ' Sign. Curv. Pred. (Blue Concave; Red Convex)'],[-1*maxcurv maxcurv]);
            colormap(mapname);
            if saveresults,G2FPredictionV2.savingResults(v,RefScan);end
            % imaging the signed curvature Differences
%             [v,RefScan] = G2FPredictionV2.cloneAndDisplay(vref,RefScanV,'Value',comparison.signedCurvatureDiff,'Indexed',[name ' Sign. Curv. Diff (Blue Decreased Convexity; Red Increased convexity)'],[]);
%             colormap(mapname);
%             if saveresults,G2FPredictionV2.savingResults(v,RefScan);end
       end
       function exportFaces(obj)
                disp('Exporting Faces..');
                exportWavefront(obj.BaseFace,'BaseFace.obj');
                exportWavefront(obj.PredFace,'PredictedFace.obj');
                exportWavefront(obj.AmplFace,'AmplifiedFace.obj');
                disp('--------------------------------------------------------------');
       end
   end
   methods % OBJECT INTERACTION
       function obj = clear(obj)
           bobj = obj;
           obj = G2FPredictionV4;
           delete(bobj);
       end
   end
   methods (Static = true)
       function [scan] = getBaseFaceStatic(A,Aval,B,Model)
            [A,B] = eliminateNAN(A,B);
            [~,~,~,~,MC] = plsregress(A,B,size(A,2));
            mC = mean(A);
            mShape = mean(B);
            DX = mC-Aval;
            DY = DX*MC(2:end,:);
            NY = mShape-DY;
            scan = getScan(Model,NY');
            scan.Material = 'Dull';scan.SingleColor = [0.8 0.8 0.8];
       end
       function [v,RefScan] = cloneAndDisplay(vref,RefScanV,what,value,color,name,range)
             if nargin < 7, range = []; end
             RefScan = clone(RefScanV);delete(RefScan.PoseLM);
             v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
             eval(['RefScan.' what ' = value;']);
             switch color
                 case 'Single'
                     RefScan.ColorMode = 'Single';
                 case 'Indexed'
                      RefScan.ColorMode = 'Indexed';
                      if isempty(range)
                        set(v.RenderAxes,'clim',[-1*max(abs(RefScan.Value)) max(abs(RefScan.Value))]);
                      else
                        set(v.RenderAxes,'clim',[range(1) range(2)]);
                      end
                      colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
                 case 'Texture'
                     RefScan.ColorMode = 'Texture';
             end
             set(v.Figure,'Name',name);set(v.Figure,'Position',get(vref.Figure,'Position'));
             v.SceneLightPosition = vref.SceneLightPosition;
             drawnow;
       end
       function savingResults(v,RefScan)
          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
       end 
   end
end


%                         if g==1 % Figure out reference to use
%                            obj.SNPRedShape = SNPCont.SNPs{g}.ST4RedShape;
%                            switch obj.SNPRedShape
%                                case 1
%                                    Ref.Coeff = BaseCont.ST4Ref.RedCoeff;
%                                    Ref.Scan = clone(BaseCont.ST4Ref.RedScan);
%                                    Ref.Cov = [];
%                                case 0
%                                    Ref.Coeff = BaseCont.ST4Ref.Coeff;
%                                    Ref.Scan = clone(BaseCont.ST4Ref.Scan);
%                                    switch BaseCont.RedShapeSpaceType
%                                         case 'X'
%                                            Ref.Cov = [BaseCont.ST4Ref.AvgCov, BaseCont.ST4Ref.AvgGB];
%                                         case 'RIP'
%                                            Ref.Cov = [BaseCont.ST4Ref.CovRIP, BaseCont.ST4Ref.GBRIP];
%                                    end
%                            end
%                            InitRef = Ref;
%                         end
