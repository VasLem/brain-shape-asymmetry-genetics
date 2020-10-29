classdef scanPCA < myPCA
    % This is an intermediate abstract class between PCA and every other
    % type of PCA class dealing with objects represented by meshes in particular faces
    % This class is storing methods that are applicable to the relevant
    % subclasses
    methods % Constructor
        function obj = scanPCA(varargin)
            obj = obj@myPCA(varargin{:});         
        end
    end
    methods % Specific Interface Functions
        function animatePC(obj,pc,range,rec)
            % animatePC(obj,pc,range,rec)
            % animate principal component, obj, pc = principal component
            % range [-x x], record true or false;
            if nargin == 1;
               pc = 1;range = [-3 3];rec = false;
            elseif nargin == 2
               range = [-3 3];rec = false;
            elseif nargin == 3
               rec = false;
            end
            MovieFile = [];
            if rec % initialize recording
               [filename, pathname, filterindex] = uiputfile( ...
                                             {'*.avi','Avi files (*.avi)'}, ...
                                              'Save as', 'Movie'); %#ok<NASGU>
               if filename == 0, 
                  rec = false; 
               else
                  cd(pathname);
                  MovieFile = avifile(filename);
                  %MovieFile = VideoWriter(filename);
                  MovieFile.fps = 3;
                  MovieFile.compression = 'none';
               end
            end
            coeff = zeros(1,obj.nrEV);
            coeff(pc) = range(1)*obj.EigStd(pc);
            [f,scan] = initializeShowPC(obj,coeff);
            in = input('Ready? y/n:','s');
            while strcmp(in,'y')
                  MovieFile = runPC(obj,pc,range,scan,rec,f,MovieFile);
                  in = input('Again? y/n:','s');
            end
            if rec, MovieFile = close(MovieFile);end
            if superClass.isH(scan),delete(scan);end
            try delete(f);catch end %#ok<CTCH>
        end
        
        function animatePCs(obj,pclist,range,rec)
            % animatePC(obj,pc,range,rec)
            % animate principal component, obj, pc = principal component
            % range [-x x], record true or false;
            if nargin == 1
               pclist = 1;range = [-3 3];rec = false;
            elseif nargin == 2
               range = [-3 3];rec = false;
            elseif nargin == 3
               rec = false;
            end
            MovieFile = [];
            if rec % initialize recording
               [filename, pathname, filterindex] = uiputfile( ...
                                             {'*.avi','Avi files (*.avi)'}, ...
                                              'Save as', 'Movie'); %#ok<NASGU>
               if filename == 0, 
                  rec = false; 
               else
                  cd(pathname);
                  MovieFile = avifile(filename);
                  %MovieFile = VideoWriter(filename);
                  MovieFile.fps = 3;
                  MovieFile.compression = 'none';
               end
            end
            nPC = length(pclist);
            for i=1:1:nPC
                disp(num2str(i));
                pc = pclist(i);
                coeff = zeros(1,obj.nrEV);
                coeff(pc) = range(1)*obj.EigStd(pc);
                if i==1, [f,scan] = initializeShowPC(obj,coeff); pause;end
                MovieFile = runPC(obj,pc,range,scan,rec,f,MovieFile);
                MovieFile = runPC(obj,pc,range,scan,rec,f,MovieFile);
            end
            if rec, MovieFile = close(MovieFile);end
            if superClass.isH(scan),delete(scan);end
            try delete(f);catch end %#ok<CTCH>
        end
        
        
        function movie = runPC(obj,pc,range,scan,rec,f,movie)
             for i=range(1):0.5:range(2)
                 coeff = zeros(1,obj.nrEV);
                 coeff(pc) = i*obj.EigStd(pc);
                 updateShowPC(obj,scan,coeff);
                 drawnow;pause(0.01);
                 if rec,movie = recordFrame(movie,f);end
             end
             for i=1:1:2
                 if rec,movie = recordFrame(movie,f);end
                 pause(0.01);
             end
             for i=range(2):-0.5:range(1)
                 coeff = zeros(1,obj.nrEV);
                 coeff(pc) = i*obj.EigStd(pc);
                 updateShowPC(obj,scan,coeff);
                 drawnow;pause(0.01);
                 if rec,movie = recordFrame(movie,f);end
             end
             function movie = recordFrame(movie,f)
                      switch class(f)
                          case 'double'
                             set(f,'Resize','off');
                             im = getframe(f);
                          case 'viewer3DObj'
                             set(f.Figure,'Resize','off');
                             im = getframe(f.Figure);
                      end                      
                      movie = addframe(movie,im);
             end
        end
        function checkAverage(obj)
            if isempty(obj.Average),obj.Average = meshObj;end
        end
        function out = getScan(obj,in)
            % function to convert any input to a scan object;
            out = embedInScan(obj,getStruc(obj,in));          
        end
    end
    methods (Abstract = true)
        [f,scan] = initializeShowPC(obj,coeff);
        updateShowPC(obj,scan,coeff);
        objout = embedInScan(obj,out);
    end
end