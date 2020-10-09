classdef batch < superClass
    % This is the superclass for batch processing and collecting
    % it is able to look for all similar type files within a folder and its
    % subfolders
    properties
        InputFolder = [];
        IncludeSubFolders = false;
        Filter = '*.*';
        gui = [];
        Interupt = false;
        StatusBar = [];
    end
    properties (Dependent = true)
        InputFiles;
        Folders;
        nrInputFiles;
        nrFolders;
    end
    methods %Constructor
        function obj = batch(varargin)
          obj = obj@superClass(varargin{:});
        end
    end
    methods% special Setting and Getting
        function out = get.nrInputFiles(obj)
            out = length(obj.InputFiles);
        end
        function out = get.nrFolders(obj)
            out = length(obj.Folders);
        end
        function out = get.Folders(obj)
            if isempty(obj.InputFolder), out = []; return; end
            out{1} = obj.InputFolder;
            if ~obj.IncludeSubFolders, return; end
            out = [out batch.getSubFolders(out{1})];
        end
        function out = get.InputFiles(obj)
            Folders = obj.Folders;
            out = {};
            if isempty(Folders), return; end
            for f=1:1:length(Folders)
                cd(Folders{f});
                tmp = dir(obj.Filter);
                if isempty(tmp), continue; end
                files = cell(1,length(tmp));
                for i=1:1:length(files)
                    files{i} = [Folders{f} '\' tmp(i).name];
                end
                out = [out files]; %#ok<AGROW>
            end  
        end
    end
    
    methods% Interface Functions
        function setInputFolder(obj)
            path = uigetdir(obj.InputFolder,'Input Folder');
            if path==0, return; end
            obj.InputFolder = path;
        end
        function killStatusBar(obj)
            try
                delete(obj.StatusBar);
            catch
            end
        end
    end
    methods (Static = true)
        function [filename,path,extension] = filestr2parts(str)
            index = strfind(str,'\');
            path = str(1:index(end));
            filename = str(index(end)+1:end);
            index = strfind(filename,'.');
            if isempty(index), extension = []; return; end
            extension = filename(index+1:end);
            filename = filename(1:index-1);
        end
        function out = FileExist(file)
                 [filename,path,extension] = batch.filestr2parts(file);
                 oldpath = pwd;
                 cd(path);
                 out = false;
                 % normal test
                 if isempty(extension)
                    tmp = dir([filename '.*']);
                 else
                    tmp = dir([filename '.' extension]);
                 end
%                  if isempty(tmp), out = false; end
%                  cd(oldpath);
                 if ~isempty(tmp), out = true;cd(oldpath);return; end
                 % cleaned up string test
                 filename = cleanUpString(filename);
                 if isempty(extension)
                    tmp = dir([filename '.*']);
                 else
                    tmp = dir([filename '.' extension]);
                 end
                 if ~isempty(tmp), out = true;cd(oldpath);end
        end
        function str = parts2filestr(filename,path,extension)
            str = [path filename '.' extension];
        end
        function subFolders = getSubFolders(Folder)
            cd(Folder);
            all = dir;
            subFolders = {};
            for i=1:1:length(all);
                %cd(Folder);
                if all(i).isdir&&~strcmp(all(i).name,'.')&&~strcmp(all(i).name,'..')
                   newFolder = [Folder '\' all(i).name];
                   subFolders = [subFolders newFolder batch.getSubFolders(newFolder)]; %#ok<AGROW>
                end
            end         
        end
    end
end