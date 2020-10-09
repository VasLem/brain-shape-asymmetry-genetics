classdef batchCollector < batchLoader
    % This is the superclass for batch processing and collecting
    % it is able to look for all similar type files within a folder and its
    % subfolders
    properties
        Data = [];
        Shape = true;
        TextureColor = false;
        TextureMap = false;
        RefScan = [];
        MaskIndex = [];
    end
    methods %Constructor
        function obj = batchCollector(varargin)
          obj = obj@batchLoader(varargin{:});
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
            disp('Starting to collect');           
        end
        function scan = innerFunction(obj,scan)
            warning('off','all');
            if ~isempty(obj.MaskIndex), crop(scan,'VertexIndex', obj.MaskIndex); end
            if obj.Shape, alignPose(scan,obj.RefScan); obj.Data.Shape(:,obj.Current) = scan.Vertices(:); end
            if obj.TextureColor, obj.Data.TextureColor(:,obj.Current) = uint8(scan.TextureColor(:)*255); end
            if obj.TextureMap, obj.Data.TextureMap(:,obj.Current) = scan.TextureMap.uint8Image(:); end
            obj.Data.Names{obj.Current} = scan.Tag;
            warning('on','all');
        end
        function done(obj)
%            cd(obj.InputFolder);
%            cd ..;
%            if ~isempty(obj.OutputFolder), cd(obj.OutputFolder); end
%            str = [obj.InputFolder '/CollectedData'];
%            if obj.Shape, str = [str '_Shape']; end
%            if obj.TextureColor, str = [str '_TextureColor']; end
%            if obj.TextureMap, str = [str '_TextureMap']; end
%            save(str,'obj');
%            clear Data;
%            disp('Done'); 
        end
    end
end