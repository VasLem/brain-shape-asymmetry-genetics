classdef G2FPredictionV2 < superClass
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
      BSNP = 2;
      R2 = 'CO';%LM CO NO
   end
   properties (Dependent = true)
      nSNPUsed; 
   end
   methods % Constructor
        function obj = G2FPredictionV2(varargin)
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
                obj.BaseFace = G2FPredictionV2.getBaseFaceStatic(BaseFaceModel.TRX,XVal,BaseFaceModel.TRBT,BaseFaceModel.Texture);
                tmp = G2FPredictionV2.getBaseFaceStatic(BaseFaceModel.TRX,XVal,BaseFaceModel.TRBS,BaseFaceModel.Shape);
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
                        case 'NO'
                            weightsR(:,g) = ones(size(weightsR(:,g)));
                        otherwise
                            delete(f);
                            error('Wrong option for obj.R2: options are LM CO or NO');
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
   methods % Visualisation functions
       function imageFacePrediction(obj,RefScanV,vref,name,saveresults)
            BaseFace = obj.BaseFace; %#ok<*PROP>
            AmplFace = obj.AmplFace;
            comparison = obj.IdAss;
            BaseFace.Material = 'Dull';
            AmplFace.Material = 'Dull';
            % imaging the scan
            [v,RefScan] = G2FPredictionV2.cloneAndDisplay(vref,RefScanV,'Vertices',BaseFace.Vertices,'Texture',[name '_Base'],[]);
            if saveresults,G2FPredictionV2.savingResults(v,RefScan);end
            % imaging the norm
            [v,RefScan] = G2FPredictionV2.cloneAndDisplay(vref,RefScanV,'Vertices',AmplFace.Vertices,'Texture',[name '_Pred'],[]);
            if saveresults,G2FPredictionV2.savingResults(v,RefScan);end
            % imaging the magnitude of the displacements
            val = vDistances(BaseFace,AmplFace);
            [v,RefScan] = G2FPredictionV2.cloneAndDisplay(vref,AmplFace,'Value',val,'Indexed',[name 'Disp (Blue NO disp; Red MAX disp)'],[0 max(val)]);
            if saveresults,G2FPredictionV2.savingResults(v,RefScan);end
            % imaging the displacements along the normals
            [v,RefScan] = G2FPredictionV2.cloneAndDisplay(vref,RefScanV,'Value',comparison.NormalDistances,'Indexed',[name ' Normal (Blue Inward disp; Red Outward disp)'],[]);
            if saveresults,G2FPredictionV2.savingResults(v,RefScan);end
            % imaging the Area Ratios
            [v,RefScan] = G2FPredictionV2.cloneAndDisplay(vref,RefScanV,'Value',comparison.AreaRatios,'Indexed',[name ' Area (Blue Decrease; Red Increase )'],[]);
            if saveresults,G2FPredictionV2.savingResults(v,RefScan);end
            % imaging the signed Curvatures
            maxcurv = min(max(abs(comparison.signedCurvature1)),max(abs(comparison.signedCurvature2)));
            [v,RefScan] = G2FPredictionV2.cloneAndDisplay(vref,BaseFace,'Value',comparison.signedCurvature2,'Indexed',[name ' Sign. Curv. Base (Blue Concave; Red Convex)'],[-1*maxcurv maxcurv]);
            if saveresults,G2FPredictionV2.savingResults(v,RefScan);end
            [v,RefScan] = G2FPredictionV2.cloneAndDisplay(vref,AmplFace,'Value',comparison.signedCurvature1,'Indexed',[name ' Sign. Curv. Pred. (Blue Concave; Red Convex)'],[-1*maxcurv maxcurv]);
            if saveresults,G2FPredictionV2.savingResults(v,RefScan);end
            % imaging the signed curvature Differences
            [v,RefScan] = G2FPredictionV2.cloneAndDisplay(vref,RefScanV,'Value',comparison.signedCurvatureDiff,'Indexed',[name ' Sign. Curv. Diff (Blue Decreased Convexity; Red Increased convexity)'],[]);
            if saveresults,G2FPredictionV2.savingResults(v,RefScan);end
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