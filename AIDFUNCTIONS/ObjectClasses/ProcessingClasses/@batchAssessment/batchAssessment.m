classdef batchAssessment < batchProcess
    properties
       Significance = 2;
       Threshold = 3;
       CompensateScale = false;
       OutliersOnly = false;
       OutliersBinary = false;
       OutliersThreshold = 0.8;
       MaskIndex = [];
       NormFolder = [];
    end
    methods %Constructor
        function obj = batchAssessment(varargin)
          obj = obj@batchProcess(varargin{:});
          obj.Overwrite = false;
          obj.Filter = '*.mat';
        end
    end
    methods% special getting and setting
    end
    methods% Interface Functions
        function done(obj)
            % dummy, do nothing, in a convertor setup
        end
        function start(obj)
            % dummy do nothing in a convertor setup
        end
        function processFunction(obj,InputFile,OutputFile) %#ok<INUSL>
            [filename,path,ext] = batch.filestr2parts(InputFile); %#ok<NASGU>
            disp(num2str(obj.Current));
            cd(obj.NormFolder);
            normfiles = dir('*.mat');
            name = cleanUpString(filename);
            %name = filename;
            hasNorm = false;
            for j=1:1:size(normfiles,1);
                normname = cleanUpString(normfiles(j).name(1:end-4));
                %normname = normfiles(j).name(1:end-4);
%                 if ~isempty(strfind(name,normname))
%                     hasNorm = true;
%                     break;
%                 end
                  if strcmp(name,normname)
                    hasNorm = true;
                    break;
                  end
            end
            if ~hasNorm,disp('NO NORM');return; end
            % Testing if the assessment already exists
            if ~obj.Overwrite
                cd(obj.OutputFolder);
                tmp = dir([name '.mat']);
                if ~isempty(tmp),return;end 
            end
            ass = assessment(obj2struc(obj));
            ass.Tag = name;
            in = load(InputFile);
            loaded = fields(in);
            ass.Scan = in.(loaded{1});
            ass.Scan.Tag = name;
            if ~isempty(obj.MaskIndex), crop(ass.Scan,'VertexIndex',obj.MaskIndex); end
            cd(obj.NormFolder)
            in = load(normfiles(j).name);
            loaded = fields(in);
            ass.Norm = in.(loaded{1});
            if ~isempty(obj.MaskIndex), crop(ass.Norm,'VertexIndex',obj.MaskIndex); end
            ass.Norm.Tag = normname;
            %update(ass);
            cd(obj.OutputFolder)
            save(ass.Tag,'ass');
        end     
    end
end