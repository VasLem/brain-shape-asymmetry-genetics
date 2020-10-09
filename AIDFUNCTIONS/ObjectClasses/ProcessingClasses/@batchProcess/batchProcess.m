classdef batchProcess < batch
    % This is the superclass for batch processing and collecting
    % it is able to look for all similar type files within a folder and its
    % subfolders
    properties
        OutputFolder = [];
        Prefix = [];
        Sufix = [];
        FailedIndex = [];
        FailedError = [];
        Overwrite = true;
        Current = 0;
        Debug = false;
    end
    properties (Dependent = true)
        Failed;
        OutputFiles;
        nrOutputFiles;
    end
    methods %Constructor
        function obj = batchProcess(varargin)
          obj = obj@batch(varargin{:});
        end
    end
    methods% special Setting and Getting
        function out = get.nrOutputFiles(obj)
            out = length(obj.InputFiles);
        end
        function out = get.OutputFiles(obj)
            InputFiles = obj.InputFiles;
            if isempty(InputFiles), out = {}; return; end
            out = InputFiles;
            for i=1:1:length(out)
                [filename,path] = batch.filestr2parts(out{i});
                if ~isempty(obj.OutputFolder), path = obj.OutputFolder; end
                out{i} = [path  obj.Prefix filename obj.Sufix];
            end
        end
        function obj = set.OutputFolder(obj,in)
            if~strcmp(in(end),'\')
              in = [in '\'];
            end
            obj.OutputFolder = in;
            if ~isdir(obj.OutputFolder),mkdir(obj.OutputFolder);end
        end
        function out = get.Failed(obj)
            out = find(obj.FailedIndex);
        end
    end  
    methods% Interface Functions
        function setOutputFolder(obj)
            path = uigetdir(obj.InputFolder,'Output Folder');
            if path==0, return; end
            obj.OutputFolder = path;
        end
        function process(obj,index)
            oldpath = pwd;
            start(obj);
            InputFiles = obj.InputFiles;
            if isempty(InputFiles),disp('Nothing to Process'), return; end
            nrfiles = length(InputFiles);
            obj.FailedIndex = zeros(1,nrfiles);
            obj.FailedError = cell(1,nrfiles);
            if nargin < 2
               index = (1:nrfiles);
            end
            %obj.StatusBar=statusbar(obj.Type);
            wb = waitbar(0,'processing');
            for i=1:1:length(index)
                obj.Current = obj.Current+1;
                if obj.Interupt, break; end
                if ~obj.Debug
                    try
                       [filename,path] = batch.filestr2parts(InputFiles{index(i)});
                       if ~isempty(obj.OutputFolder), path = obj.OutputFolder; end
                       OutputFile = [path  obj.Prefix filename obj.Sufix];
                       set(obj.StatusBar,'Name',filename);
                       drawnow;processFunction(obj,InputFiles{index(i)},OutputFile);       
                    catch the_error
                       obj.FailedIndex(index(i)) = 1;
                       obj.FailedError{index(i)} = the_error;
                    end
                else
                    [filename,path] = batch.filestr2parts(InputFiles{index(i)});
                    if ~isempty(obj.OutputFolder), path = obj.OutputFolder; end
                    OutputFile = [path  obj.Prefix filename obj.Sufix];
                    set(obj.StatusBar,'Name',filename);drawnow;
                    processFunction(obj,InputFiles{index(i)},OutputFile);
                end
                %statusbar((i/length(index)),obj.StatusBar);
                waitbar(i/length(index),wb);drawnow;
            end
            killStatusBar(obj);
            obj.Current = 0;
            cd(oldpath);
            done(obj);
            delete(wb);
        end
        function showError(obj,index)
            obj.FailedError{index}
            for i=1:1:length(obj.FailedError{index}.stack)
                obj.FailedError{index}.stack(i)
            end           
        end
    end
    methods (Abstract = true);
        start(obj);
        processFunction(obj,InputFile,OutputFile);
        done(obj);
    end
end