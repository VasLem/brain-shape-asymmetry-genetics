classdef batchPreMapper < batchConvertor
    properties       
       CenterVertices = true;
       IndicatePoseLM = true;
       AdjustScale = false;
       ScaleRefScan = [];
       Scale = false;
       ScaleFactor = 1;
       DeleteBorder = false;
       DeleteBorderTimes = 10;
       Crop = false;
       CropDistance = 200;
       PreView = false;
       PostView = false;
       RbfFit = false;
       DownSampleTextureMap = false;
       DownSampleFaces = false;
       DownSampleFacesPercentage = 60;
       CleanManifold = false;
    end
    methods %Constructor
        function obj = batchPreMapper(varargin)
          obj = obj@batchConvertor(varargin{:});
          obj.InputFormat = 'Mat Object';
          obj.OutputFormat = 'Mat Object';
          obj.Overwrite = false;
        end
    end
    methods% special getting and setting
        function out = get.ScaleRefScan(obj)
            out = obj.ScaleRefScan;
            if ~superClass.isH(out), out = []; end
        end
    end
    methods% Interface Functions
        function scan = innerFunction(obj,scan)
            %if ~isempty(scan.PoseLM), return; end
            if obj.CenterVertices, centerVertices(scan); end
            if obj.Scale
                 T = scaledRigidTM;
                 T.Scale = obj.ScaleFactor;
                 transform(T,scan);
            end
            if obj.CleanManifold, cleanManifold(scan); end
            if obj.PreView, batchView(scan); end
            if obj.DeleteBorder, deleteBorder(scan,obj.DeleteBorderTimes); end
            if obj.DownSampleFaces, reduceTriangles(scan,obj.DownSampleFacesPercentage/100); end
            if obj.IndicatePoseLM, indicatePoseLM(scan); end
            if obj.AdjustScale
               T = scaledRigidTM;
               if isempty(obj.ScaleRefScan),
                  scale = 1;
               else
                  match(T,scan.PoseLM,obj.ScaleRefScan.PoseLM);
                  scale = max(1,round(T.Scale));
               end
               T.Scale = 1/scale;
               transform(T,scan);
            end
            if obj.Crop, poseLMCrop(scan,obj.CropDistance);end
            if obj.PostView, batchView(scan); end
            if obj.DownSampleTextureMap
                try
                 factor = round(numel(scan.TextureMap.Image)/4000000);
                 if factor> 0, downSample(scan.TextureMap,factor); end
                catch
                end
            end
            if obj.CenterVertices, centerVertices(scan); end
            %if obj.RbfFit, rbfFit(scan); end
            % dummy, do nothing
        end     
    end
end

% classdef batchPreMapper < batchConvertor
%     properties       
%        CenterVertices = true;
%        IndicatePoseLM = true;
%        AdjustScale = false;
%        ScaleRefScan = [];
%        DeleteBorder = false;
%        DeleteBorderTimes = 10;
%        Crop = false;
%        CropDistance = 200;
%        PreView = false;
%        PostView = false;
%        RbfFit = false;
%        DownSampleTextureMap = false;
%        DownSampleFaces = false;
%        DownSampleFacesPercentage = 60;
%     end
%     methods %Constructor
%         function obj = batchPreMapper(varargin)
%           obj = obj@batchConvertor(varargin{:});
%           obj.InputFormat = 'Mat Object';
%           obj.OutputFormat = 'Mat Object';
%           obj.Overwrite = false;
%         end
%     end
%     methods% special getting and setting
%         function out = get.ScaleRefScan(obj)
%             out = obj.ScaleRefScan;
%             if ~superClass.isH(out), out = []; end
%         end
%     end
%     methods% Interface Functions
%         function scan = innerFunction(obj,scan)
%             %if ~isempty(scan.PoseLM), return; end
%             if obj.CenterVertices, centerVertices(scan); end
%             if obj.PreView, batchView(scan); end
%             if obj.DeleteBorder, deleteBorder(scan,obj.DeleteBorderTimes); end
%             if obj.DownSampleFaces, reduceTriangles(scan,obj.DownSampleFacesPercentage/100); end
%             if obj.IndicatePoseLM, indicatePoseLM(scan); end
%             if obj.AdjustScale
%                T = scaledRigidTM;
%                if isempty(obj.ScaleRefScan),
%                   scale = 1;
%                else
%                   match(T,scan.PoseLM,obj.ScaleRefScan.PoseLM);
%                   scale = max(1,round(T.Scale));
%                end
%                T.Scale = 1/scale;
%                transform(T,scan);
%             end
%             %if obj.Crop, poseLMCrop(scan,obj.CropDistance);end
%             if obj.PostView, batchView(scan); end
%             if obj.DownSampleTextureMap
%                 try
%                  factor = round(numel(scan.TextureMap.Image)/4000000);
%                  if factor> 0, downSample(scan.TextureMap,factor); end
%                 catch
%                 end
%             end
%             %if obj.RbfFit, rbfFit(scan); end
%             % dummy, do nothing
%         end     
%     end
% end