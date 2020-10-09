classdef G2FPredictionV4 < superClass
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
      SNP;
      SNPUsed = [];
      RIPG;
      WithRIPG = false;
      BaseFace;
      PredFace;
      AmplFace;
      AF = 4;
      IdAss;
      BSNP = 1;
      R2 = 'LM';
   end
   properties (Dependent = true)
      nSNPUsed; 
   end
   methods % Constructor
        function obj = G2FPredictionV4(varargin)
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
           if datenum(date)>datenum('31-Oct-2015'), error('License Expired, please contact peter.claes@esat.kuleuven.be');end
           if datenum(date)+30>datenum('31-Oct-2015')
              warning(['License will expire in ' num2str(datenum('31-Oct-2015')-datenum(date)) ' days , please contact peter.claes@esat.kuleuven.be']);
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
           obj.AncNames = raw(9:8+obj.nAnc,1);
           obj.SNPNames = raw(9+obj.nAnc:8+obj.nAnc+obj.nSNP,1);
           tmp = raw(9+obj.nAnc:8+obj.nAnc+obj.nSNP,2);
           obj.SNP = nan(size(tmp));
           for i=1:1:length(tmp)
               if isnumeric(tmp{i})
                  obj.SNP(i) = tmp{i};
               end
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
                XVal = [obj.Sex obj.Age obj.BMI obj.H obj.Anc'];
                index = find(isnan(XVal));
                if ~isempty(index), XVal(index) = XAvg(index);end
                obj.BaseFace = G2FPredictionV3.getBaseFaceStatic(BaseFaceModel.TRX,XVal,BaseFaceModel.TRBT,BaseFaceModel.Texture);
                tmp = G2FPredictionV3.getBaseFaceStatic(BaseFaceModel.TRX,XVal,BaseFaceModel.TRBS,BaseFaceModel.Shape);
                obj.BaseFace.Vertices = tmp.Vertices;
                obj.BaseFace.ColorMode = 'Texture';
                disp('--------------------------------------------------------------');
       end
       function overlaySNPs(obj,SNPs)
                if isempty(obj.BaseFace), error('First create BaseFace or call getFacialPrediction');end
                disp('Matching Genotypes & overlaying effects...');
                f = statusbar('Processing...');drawnow;
                good = zeros(1,obj.nSNP);
                shapeR = nan*zeros(obj.BaseFace.nrV*3,obj.nSNP);
                weightsR = nan*zeros(obj.BaseFace.nrV*3,obj.nSNP);
                for g=1:1:obj.nSNP
                    if isnan(obj.SNP(g)), disp('SNP = nan'); continue; end
                    i = find(strcmp(obj.SNPNames{g},SNPs.Names));
                    if isempty(i), disp('SNP not found'); continue; end
                    if ~SNPs.Succes; continue; end
                    if SNPs.TestInd(i)==0, disp('Bad SNP Test');continue; end
                    switch obj.SNP(g)
                        case -1
                            if isempty(SNPs.AAMorphs{i}), disp('Empty AA'); continue;end
                            switch obj.BSNP
                                case 1
                                    morph = clone(SNPs.AAMorphs{i});
                                case 2
                                    morph = clone(SNPs.AAMorphs2{i});
                                otherwise
                                    morph = clone(SNPs.AAMorphs3{i});
                            end
                        case 0
                            if isempty(SNPs.ABMorphs{i}), disp('Empty AB'); continue;end
                            switch obj.BSNP
                                case 1
                                    morph = clone(SNPs.ABMorphs{i});
                                case 2
                                    morph = clone(SNPs.ABMorphs2{i});
                                otherwise
                                    morph = clone(SNPs.ABMorphs3{i});
                            end
                        case 1
                            if isempty(SNPs.BBMorphs{i}), disp('Empty BB'); continue;end
                            switch obj.BSNP
                                case 1
                                    morph = clone(SNPs.BBMorphs{i});
                                case 2
                                    morph = clone(SNPs.BBMorphs2{i});
                                otherwise
                                    morph = clone(SNPs.BBMorphs3{i});
                            end
                        otherwise
                            % unknown genotype
                            continue;
                    end
                    diff = morph.Vertices-SNPs.RefScan{i}.Vertices+obj.BaseFace.Vertices;
                    shapeR(:,g) = diff(:);
                    switch obj.R2
                        case 'LM'
                            tmp = repmat(SNPs.R2(i,:),3,1);
                            weightsR(:,g) = tmp(:);
                        case 'CO'
                            weightsR(:,g) = SNPs.R2CO(i,:)';
                    end
                    good(g) = 1;
                    statusbar(g/obj.nSNP,f);drawnow;
                end
                shapeR = shapeR(:,good==1);
                weightsR = weightsR(:,good==1);
                shaperes = sum(shapeR.*weightsR,2)./sum(weightsR,2);
                obj.PredFace = clone(obj.BaseFace);
                obj.PredFace.Vertices = reshape(shaperes,3,obj.BaseFace.nrV);
                obj.PredFace.Value = vDistances(obj.PredFace,obj.BaseFace);
                obj.SNPUsed = find(good);
                disp([num2str(obj.nSNPUsed) ' SNPs Used']);
                delete(f);
                disp('--------------------------------------------------------------');
       end
       function getFacialPrediction(obj,BaseFaceModel,SNPs)
                getBaseFace(obj,BaseFaceModel);
                overlaySNPs(obj,SNPs);
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
   end
end