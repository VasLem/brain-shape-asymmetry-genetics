classdef batchAssessmentAverager < batchProcess
    % This is the superclass for batch processing and collecting
    % it is able to look for all similar type files within a folder and its
    % subfolders
    properties
        Average = [];
        Counter = 0;
    end
    methods %Constructor
        function obj = batchAssessmentAverager(varargin)
          obj = obj@batchProcess(varargin{:});
          obj.Filter = '*.mat';
        end
    end
    methods% special Setting and Getting
        function out = get.Average(obj)
            out = obj.Average;
            if ~superClass.isH(out), out = []; end
        end
    end  
    methods% Interface Functions
        function start(obj)
            disp('Initializing');
            obj.Average = [];
            obj.Counter = 0;
            disp('Starting to Average');           
        end
       function processFunction(obj,InputFile,OutputFile) %#ok<INUSL>
            in = load(InputFile);
            loaded = fields(in);
            ass = in.(loaded{1});
            if ass.Update, return; end
            if obj.Counter == 0
                obj.Average = clone(ass);
                if~isempty(obj.Average.Scan.TextureMap), delete(obj.Average.Scan.TextureMap); end
                if~isempty(obj.Average.Scan.PoseLM), delete(obj.Average.Scan.PoseLM); end
                if~isempty(obj.Average.Scan.Border), delete(obj.Average.Scan.Border); end
                obj.Average.Scan.TextureColor = [];obj.Average.Scan.UV = [];
                obj.Average.Scan.ColorMode = 'Single';
                if~isempty(obj.Average.Norm.TextureMap), delete(obj.Average.Norm.TextureMap); end
                if~isempty(obj.Average.Norm.PoseLM), delete(obj.Average.Norm.PoseLM); end
                if~isempty(obj.Average.Norm.Border), delete(obj.Average.Norm.Border); end
                obj.Average.Norm.TextureColor = [];obj.Average.Norm.UV = [];
                obj.Average.Norm.ColorMode = 'Single';
            else
                out =  alignPose(ass.Scan,obj.Average.Scan,nan); 
                obj.Average.Scan.Vertices = obj.Average.Scan.Vertices+out.Vertices;
                delete(out);
                out =  alignPose(ass.Norm,obj.Average.Norm,nan); 
                obj.Average.Norm.Vertices = obj.Average.Norm.Vertices+out.Vertices;
                delete(out);
                obj.Average.Dysmorphogram = obj.Average.Dysmorphogram+ass.Dysmorphogram;
                obj.Average.PercDysmorph = obj.Average.PercDysmorph+ass.PercDysmorph;
                obj.Average.ThresholdMap = obj.Average.ThresholdMap+ass.ThresholdMap;
                obj.Average.PercThresh = obj.Average.PercThresh+ass.PercThresh;
                obj.Average.DistanceMap = obj.Average.DistanceMap+ass.DistanceMap;
                obj.Average.RMSE = obj.Average.RMSE+ass.RMSE;
            end
            obj.Counter = obj.Counter+1;
        end
        function done(obj)
         disp('Done');
         if isempty(obj.Average), return; end
         obj.Average.Scan.Vertices = obj.Average.Scan.Vertices/obj.Counter;
         obj.Average.Norm.Vertices = obj.Average.Norm.Vertices/obj.Counter;
         obj.Average.Dysmorphogram = obj.Average.Dysmorphogram/obj.Counter;
         obj.Average.PercDysmorph = obj.Average.PercDysmorph/obj.Counter;
         obj.Average.ThresholdMap = obj.Average.ThresholdMap/obj.Counter;
         obj.Average.PercThresh = obj.Average.PercThresh/obj.Counter;
         obj.Average.DistanceMap = obj.Average.DistanceMap/obj.Counter;
         obj.Average.RMSE = obj.Average.RMSE/obj.Counter;
         obj.Average.Tag = 'Average';
         obj.Average.Scan.Tag = 'Average Scan';
         obj.Average.Norm.Tag = 'Average Norm';
         obj.Average.Norm.Selected = true;
         centerVertices(obj.Average.Norm);
         centerVertices(obj.Average.Scan);
        end
    end
end