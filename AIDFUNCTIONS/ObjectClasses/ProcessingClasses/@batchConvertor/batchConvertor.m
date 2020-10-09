classdef batchConvertor < batchLoader
    % This is the superclass for batch processing and collecting
    % it is able to look for all similar type files within a folder and its
    % subfolders
    properties
        OutputFormat = [];
        ExportFormats = {'Mat Object','Mat Struct','Obj Wavefront','Off','Mat Custom'};
    end
    properties (Dependent = true)
        OutputExtension;
    end
    methods %Constructor
        function obj = batchConvertor(varargin)
          obj = obj@batchLoader(varargin{:});
          obj.InputFormat = 'Obj 3DMD(2Pod)';
          obj.OutputFormat = 'Mat Object';
        end
    end
    methods% special Setting and Getting
        function obj = set.OutputFormat(obj,in)
            switch in
                case obj.ExportFormats
                otherwise
                    error('Unknown output format');
            end
            obj.OutputFormat = in;
        end
        function out = get.OutputExtension(obj)
            switch obj.OutputFormat
                case {'Mat Object' 'Mat Struct' 'FMM Original' 'FMM Mapped' 'Mat Custom'}
                    out = '.mat';
                case {'Obj Wavefront' 'Obj 3DMD(2Pod)'}
                    out = '.obj';
                case 'Wrl'
                    out = '.wrl';
                case 'Smf'
                    out = '.smf';
                case 'Ply'
                    out = '.ply';
                case 'Off'
                    out = '.off';
                otherwise
                    error('Unknown input format');
            end      
        end
    end  
    methods% Interface Functions
        % make this an abstract function
        function processFunction(obj,InputFile,OutputFile)
            OutputFile = [OutputFile obj.OutputExtension];
            if ~obj.Overwrite&&batch.FileExist(OutputFile), return; end
            scan = importFunction(obj,InputFile);
            scan = innerFunction(obj,scan);
            exportFunction(obj,scan,OutputFile);
            try
                delete(scan);
            catch
                clear scan;
            end
        end
        function scan = innerFunction(obj,scan)
            % dummy, do nothing, in a convertor setup
        end
        function done(obj)
            % dummy, do nothing, in a convertor setup
        end
        function start(obj)
            % dummy do nothing in a convertor setup
        end
        function exportFunction(obj,scan,OutputFile)
            [filename,path] = batch.filestr2parts(OutputFile);
            cd(path);
            switch obj.OutputFormat
                case 'Mat Object'
                    exportObject(scan,filename);
                case 'Mat Struct'
                    exportStruct(scan,filename);
                case 'Obj Wavefront'
                    exportWavefront(scan,filename);
                case 'Off'
                    exportOff(scan,filename);
                case 'Mat Custom'
                    save(filename,'scan');
                otherwise
                    error('unknown export format');     
            end
        end
    end
end