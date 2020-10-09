classdef batchLoader < batchProcess
    % This is the superclass for batch processing and collecting
    % it is able to look for all similar type files within a folder and its
    % subfolders
    properties
        InputFormat = [];
        ImportFormats = {'Mat Object','Mat Struct','Obj Wavefront','Obj 3DMD(2Pod)','Obj 3DMD(5Pod)','Obj MevisLab','Obj Eyetronics','Wrl','Smf','Ply','Off','Nrf','FMM Original','FMM Mapped','Mat Assessment'};
    end
    properties (Dependent = true)
    end
    methods %Constructor
        function obj = batchLoader(varargin)
          obj = obj@batchProcess(varargin{:});
          obj.InputFormat = 'Mat Object';
        end
    end
    methods% special Setting and Getting
        function obj = set.InputFormat(obj,in)
            switch in
                case {'Mat Object' 'Mat Struct' 'FMM Original' 'FMM Mapped' 'Mat Assessment'}
                    obj.Filter = '*.mat';
                case {'Obj Wavefront' 'Obj 3DMD(2Pod)' 'Obj 3DMD(5Pod)' 'Obj MevisLab' 'Obj Eyetronics'}
                    obj.Filter = '*.obj';
                case 'Wrl'
                    obj.Filter = '*.wrl';
                case 'Smf'
                    obj.Filter = '*.smf';
                case 'Ply'
                    obj.Filter = '*.ply';
                case 'Off'
                    obj.Filter = '*.off';
                case 'Nrf'
                    obj.Filter = '*.nrf';
                otherwise
                    error('Unknown Input format');
            end
            obj.InputFormat = in;
        end
    end  
    methods% Interface Functions
        % make this an abstract function
        function processFunction(obj,InputFile,OutputFile)
            scan = importFunction(obj,InputFile);
            scan = innerFunction(obj,scan);
            delete(scan);
        end
        function scan = importFunction(obj,InputFile)
            [filename,path,ext] = batch.filestr2parts(InputFile);  
            cd(path);
            scan = meshObj;
            switch obj.InputFormat
                case 'Mat Object'
                    importObject(scan,[filename '.' ext]);
                    %scan.Tag = filename;
                case 'Mat Struct'
                    importStruct(scan,[filename '.' ext]);
                case 'Obj Wavefront'
                    importWavefront(scan,[filename '.' ext],'Type','Wavefront');
                case 'Obj 3DMD(2Pod)'
                    importWavefront(scan,[filename '.' ext],'Type','3DMD 2pod');
                case 'Obj 3DMD(5Pod)'
                    importWavefront(scan,[filename '.' ext],'Type','3DMD 5pod');
                case 'Obj Eyetronics'
                    importWavefront(scan,[filename '.' ext],'Type','Eyetronics');
                case 'Obj MevisLab'
                    importWavefrontMevisLab(scan,[filename '.' ext],'Type','MevisLab');    
                case 'Wrl'
                    importWrl(scan,[filename '.' ext]);
                case 'Smf'
                    importSmf(scan,[filename '.' ext]);
                case 'Ply'
                    importPly(scan,[filename '.' ext]);
                case 'Off'
                    importOff(scan,[filename '.' ext]);
                case 'Nrf'
                    importNrf(scan,[filename '.' ext]);
                case 'FMM Original'
                    importFMM(scan,[filename '.' ext],'Type','FMM Original');
                case 'FMM Mapped'
                    importFMM(scan,[filename '.' ext],'Type','FMM Mapped');
                case 'Mat Assessment'
                    delete(scan);
                    in = load(filename);
                    loaded = fields(in);
                    scan = assessment(in.(loaded{1}));
                otherwise
                    error('unknown input format');
            end
            scan.Tag = filename;
        end
    end
    methods (Abstract = true)
            scan = innerFunction(obj,scan);    
    end
end