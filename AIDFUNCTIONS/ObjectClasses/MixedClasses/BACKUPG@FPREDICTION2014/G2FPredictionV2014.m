classdef G2FPredictionV2014 < superClass
   properties
      ID;
      Sex = nan;
      Age = nan;
      BMI= nan;
      H = nan;
      nAnc;
      AncNames;
      Anc
      nSNP
      SNPNames;
      SNPG;
      SNPUsed = [];
      %nSNPUsed;
      RIPG;
      WithRIPG = false;
      %RIPS;
      %RIPSE;
      %RIPA;
      %RIPAE;
      BaseFace;
      PredFace;
      PredFaceOld;
      AmplFace;
      AF = 4;
      IdAss;
   end
   properties (Dependent = true)
      nSNPUsed; 
   end
   methods % Constructor
        function obj = G2FPredictionV2014(varargin)
            obj = obj@superClass(varargin{:});         
        end
   end
   methods % special getting and setting
        function out = get.nSNPUsed(obj)
           out = length(obj.SNPUsed);
        end
   end
   methods % Interface functions
       function parsePredictionInput(obj,input)
           disp('Parsing Prediction Input...');
           if datenum(date)>datenum('31-May-2015'), error('License Expired, please contact peter.claes@esat.kuleuven.be');end
           if datenum(date)+30>datenum('31-May-2015')
              warning(['License will expire in ' num2str(datenum('31-May-2015')-datenum(date)) ' days , please contact peter.claes@esat.kuleuven.be']);
           end
           [~,~,raw] = xlsread(input);
           obj.nAnc = raw{1,2};
           obj.nSNP = raw{2,2};
           obj.WithRIPG = logical(raw{3,2});
           obj.ID = raw{4,2};
           obj.Sex = raw{5,2};
           obj.Age = raw{6,2};
           obj.H = raw{7,2};
           obj.BMI = raw{8,2};
           tmp = raw(9:8+obj.nAnc,2);
           obj.Anc = nan(size(tmp));
           for i=1:1:length(tmp)
               if isnumeric(tmp{i})
                  obj.Anc(i) = tmp{i};
               end
           end    
           %obj.Anc = cell2mat(raw(9:8+obj.nAnc,2));
           obj.AncNames = raw(9:8+obj.nAnc,1);
           obj.SNPNames = raw(9+obj.nAnc:8+obj.nAnc+obj.nSNP,1);
           tmp = raw(9+obj.nAnc:8+obj.nAnc+obj.nSNP,2);
           obj.SNPG = nan(size(tmp));
           for i=1:1:length(tmp)
               if isnumeric(tmp{i})
                  obj.SNPG(i) = tmp{i};
               end
           end
           %obj.SNPG = cell2mat(raw(9+obj.nAnc:8+obj.nAnc+obj.nSNP,2));
           if obj.WithRIPG
               tmp = raw(9+obj.nAnc:8+obj.nAnc+obj.nSNP,3);
               obj.RIPG = nan(size(tmp));
               for i=1:1:length(tmp)
                   if isnumeric(tmp{i})
                      obj.RIPG(i) = tmp{i};
                   end
               end   
               %obj.RIPG = cell2mat(raw(9+obj.nAnc:8+obj.nAnc+obj.nSNP,3));
           end
           disp('--------------------------------------------------------------');
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
       function getBaseFace(obj,BaseFaceModel)
                disp('Creating Base Face...');
                XAvg = nanmean(BaseFaceModel.TRX);
                XVal = [obj.Sex obj.Age obj.BMI obj.H obj.Anc(BaseFaceModel.Axindex)'];
                index = find(isnan(XVal));
                if ~isempty(index), XVal(index) = XAvg(index);end
                obj.BaseFace = getBaseFace(BaseFaceModel.TRX,XVal,BaseFaceModel.TRB,BaseFaceModel.AppModel);
                disp('--------------------------------------------------------------');
       end
       function overlaySNPsOld(obj,SNPs)
                if isempty(obj.BaseFace), error('First create BaseFace or call getFacialPrediction');end
                disp('Matching Genotypes & overlaying effects...');
                f = statusbar('Processing...');drawnow;
                good = [];
                shapeR = nan*zeros(obj.BaseFace.nrV*3,obj.nSNP);
                weightsR = nan*zeros(obj.BaseFace.nrV*3,obj.nSNP);
                for g=1:1:obj.nSNP
                    if isnan(obj.SNPG(g)), continue; end
                    if obj.SNPG(g)==0, continue; end
                    i = find(strcmp(obj.SNPNames{g},SNPs.Names));
                    if isempty(i), continue; end
                    switch obj.SNPG(g)
                        case -1
                            morph = SNPs.NegMorphs(:,:,i);
                        case 1
                            morph = SNPs.PosMorphs(:,:,i);
                    end
                    morph = morph - SNPs.RefFace.Vertices + obj.BaseFace.Vertices;
                    shapeR(:,g) = morph(:);
                    tmp = repmat(SNPs.LocalS(:,i)',3,1);
                    weightsR(:,g) = tmp(:);
                    good = [good g];
                    statusbar(g/obj.nSNP,f);drawnow;
                end
                shapeR = shapeR(:,good);
                weightsR = weightsR(:,good);
                shaperes = sum(shapeR.*weightsR,2)./sum(weightsR,2);
                obj.PredFaceOld = clone(obj.BaseFace);
                obj.PredFaceOld.Vertices = reshape(shaperes,3,obj.BaseFace.nrV);
                obj.PredFaceOld.Value = vDistances(obj.PredFaceOld,obj.BaseFace);
                obj.SNPUsed = good;
                disp([num2str(obj.nSNPUsed) ' SNPs Used']);
                delete(f);
                disp('--------------------------------------------------------------');     
       end
       function overlaySNPs(obj,SNPs)
                if isempty(obj.BaseFace), error('First create BaseFace or call getFacialPrediction');end
                disp('Matching Genotypes & overlaying effects...');
                f = statusbar('Processing...');drawnow;
                good = [];
                shapeR = nan*zeros(obj.BaseFace.nrV*3,obj.nSNP);
                weightsR = nan*zeros(obj.BaseFace.nrV*3,obj.nSNP);
                for g=1:1:obj.nSNP
                    if isnan(obj.RIPG(g)), continue; end
                    i = find(strcmp(obj.SNPNames{g},SNPs.Names));
                    if isempty(i), continue; end
                    FD = SNPs.PosMorphs(:,:,i)-SNPs.RefFace.Vertices;
                    rips = SNPs.RIPG(:,i);
                    m = mean(rips);s = std(rips);
                    RD = 2*(m+3*s)-m;
                    FD = FD/RD;
                    morph = obj.BaseFace.Vertices;
                    morph = morph+(2*obj.RIPG(g)-m)*FD;
                    shapeR(:,g) = morph(:);
                    tmp = repmat(SNPs.LocalS(:,i)',3,1);
                    weightsR(:,g) = tmp(:);
                    good = [good g];
                    statusbar(g/obj.nSNP,f);drawnow;
                end
                shapeR = shapeR(:,good);
                weightsR = weightsR(:,good);
                shaperes = sum(shapeR.*weightsR,2)./sum(weightsR,2);
                obj.PredFace = clone(obj.BaseFace);
                obj.PredFace.Vertices = reshape(shaperes,3,obj.BaseFace.nrV);
                obj.PredFace.Value = vDistances(obj.PredFace,obj.BaseFace);
                obj.SNPUsed = good;
                disp([num2str(obj.nSNPUsed) ' SNPs Used']);
                delete(f);
                disp('--------------------------------------------------------------');
       end
       function getFacialPrediction(obj,BaseFaceModel,SNPs)
                getBaseFace(obj,BaseFaceModel);
                overlaySNPs(obj,SNPs);
       end
       function getFacialPredictionOld(obj,BaseFaceModel,SNPs)
                getBaseFace(obj,BaseFaceModel);
                overlaySNPsOld(obj,SNPs);
                obj.PredFace = obj.PredFaceOld;
                obj.AF = 4;
       end
   end
end