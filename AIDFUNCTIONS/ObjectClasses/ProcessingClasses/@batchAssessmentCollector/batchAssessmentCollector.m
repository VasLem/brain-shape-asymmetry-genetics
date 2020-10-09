classdef batchAssessmentCollector < batchLoader
    % This is the superclass for batch processing and collecting
    % it is able to look for all similar type files within a folder and its
    % subfolders
    properties
        Data = [];
        Shape = true;
        NormShape = true;
        Dysmorphogram = true;
        BinaryDysmorphogram = true;
        BT = 0.25;
        DistanceMap = true;
        ThresholdMap = true;
        NoiseLevel = true;
        RefScan = [];
        TextureColor = false;
        TextureMap = false;
    end
    methods %Constructor
        function obj = batchAssessmentCollector(varargin)
          obj = obj@batchLoader(varargin{:});
          obj.InputFormat = 'Mat Assessment';
        end
    end
    methods% special Setting and Getting
        function out = get.RefScan(obj)
            out = obj.RefScan;
            if ~superClass.isH(out), out = []; end
        end
    end  
    methods% Interface Functions
        function start(obj)
            disp('Initializing');
            obj.Data = [];
            if isempty(obj.RefScan), error('First set refscan before collecting'); end
            nr = obj.nrInputFiles;
            obj.Data.Names = cell(1,nr);
            if obj.Shape,obj.Data.Shape = nan*zeros(obj.RefScan.nrV*3,nr);end
            if obj.TextureColor,obj.Data.TextureColor = uint8(zeros(obj.RefScan.nrV*3,nr));end
            if obj.TextureMap,obj.Data.TextureMap = uint8(zeros(numel(obj.RefScan.TextureMap.Image),nr));end
            if obj.NormShape,obj.Data.NormShape = nan*zeros(obj.RefScan.nrV*3,nr);end
            if obj.Dysmorphogram
               obj.Data.Dysmorphogram = nan*zeros(obj.RefScan.nrV,nr);
               obj.Data.PercDysmorph = nan*zeros(1,nr);
            end
            if obj.BinaryDysmorphogram
               obj.Data.BinaryDysmorphogram = nan*zeros(obj.RefScan.nrV,nr);
               obj.Data.BinaryPercDysmorph = nan*zeros(1,nr);
            end
            if obj.DistanceMap
                obj.Data.DistanceMap = nan*zeros(obj.RefScan.nrV,nr);
                obj.Data.RMSE = nan*zeros(1,nr);
            end
            if obj.ThresholdMap
               obj.Data.ThresholdMap = nan*zeros(obj.RefScan.nrV,nr);
               obj.Data.PercThresh = nan*zeros(1,nr);
            end
            if obj.NoiseLevel
                obj.Data.NoiseLevel = nan*zeros(1,nr);
            end
            disp('Starting to collect');           
        end
        function ass = innerFunction(obj,ass)
            warning('off','all');
            obj.Data.Names{obj.Current} = ass.Tag; 
            if obj.Shape, alignPose(ass.Scan,obj.RefScan); obj.Data.Shape(:,obj.Current) = ass.Scan.Vertices(:); end
            if obj.TextureColor
               try 
                  ass.Scan.TextureMap.PixInterp = '*Linear';
               catch
               end
               obj.Data.TextureColor(:,obj.Current) = uint8(ass.Scan.TextureColor(:)*255); 
            end
            if obj.TextureMap, obj.Data.TextureMap(:,obj.Current) = ass.Scan.TextureMap.uint8Image(:); end
            if obj.NormShape, alignPose(ass.Norm,obj.RefScan); obj.Data.NormShape(:,obj.Current) = ass.Norm.Vertices(:); end
            if obj.Dysmorphogram
               obj.Data.Dysmorphogram(:,obj.Current) = ass.Dysmorphogram(:); 
               obj.Data.PercDysmorph(obj.Current) = ass.PercDysmorph;
            end
            if obj.BinaryDysmorphogram
               BDtemp = ass.Dysmorphogram(:) ;
               BD = zeros(size(BDtemp));
               BD(find(BDtemp>=obj.BT)) = 1; %#ok<FNDSB>
               obj.Data.BinaryDysmorphogram(:,obj.Current) = BD; 
               obj.Data.BinaryPercDysmorph(obj.Current) = (sum(BD)/length(BD))*100;
            end 
            if obj.DistanceMap
               obj.Data.DistanceMap(:,obj.Current) = ass.DistanceMap(:);
               obj.Data.RMSE(obj.Current) = ass.RMSE;
            end
            if obj.ThresholdMap
               obj.Data.ThresholdMap(:,obj.Current) = ass.ThresholdMap(:);
               obj.Data.PercThresh(obj.Current) = ass.PercThresh;
            end
            if obj.NoiseLevel
                obj.Data.NoiseLevel(obj.Current) = ass.NoiseLevel;
            end
            obj.Data.Names{obj.Current} = ass.Tag;           
            warning('on','all');
        end
        function done(obj)
        end
    end
end