classdef propertyPCA < scanPCA
    properties
        Average = [];% Property values will now be stored in the average
        Model = [];% A PCA model related with scans, can be shape, texture or appearance
        PropNames = [];
    end
    properties(Dependent = true)
        AvgVec;% The average vertices in vector notation
        AvgProp;% The average properties
        RefScan;% The reference scan defining the Faces
        nrMC;% Number of Model Coefficients
        nrP;% Number of Properties
        nrC;% Total number of Coefficients 
    end
    methods % Constructor
        function obj = propertyPCA(varargin)
            obj = obj@scanPCA(varargin{:});         
        end
    end
    methods % Special Setting & Getting
        function out = get.Average(obj)
            checkModel(obj);
            out = obj.Model.Average;
        end
        function obj = set.Average(obj,in)
           checkModel(obj);
           obj.Model.Average = in;
        end
        function out = get.Model(obj)
            out = obj.Model;
            if ~superClass.isH(out), out = []; end
        end
        function obj = set.Model(obj,in)
            if ~isempty(obj.Model)&&~(in==obj.Model), delete(obj.Model); end
            obj.Model = in;
        end     
        function obj = set.RefScan(obj,in)
            checkModel(obj);obj.Model.RefScan = in;
        end
        function out = get.RefScan(obj)
            out = obj.Average;
        end
        function out = get.AvgVec(obj)
            if isempty(obj.Model), out = []; return; end
            if ~isfield(obj.Model.Average.UserData,'Prop'), out = []; return; end
            out = zeros(obj.nrMC,1);
            out = [out;obj.Average.UserData.Prop];
        end
        function obj = set.AvgVec(obj,in)
            obj.Average.UserData.Prop = in(obj.nrMC+1:end);
            % do nothing
        end
        function obj = set.AvgProp(obj,in)
            obj.Average.UserData.Prop = in;
        end
        function out = get.AvgProp(obj,in)
            if isempty(obj.Average), out = []; return;end
            out = obj.Average.UserData.Prop;
        end     
        function out = get.nrMC(obj)
            if isempty(obj.Model), out  =[]; return;end
            out = obj.Model.nrEV;
        end
        function out = get.nrP(obj)
            out = length(obj.PropNames);
        end
        function out = get.nrC(obj)
            out = obj.nrMC+obj.nrP;
        end
    end
    methods % Specific Interface Functions
        function checkModel(obj)
            if isempty(obj.Model), obj.Model = shapePCA; end
        end
        function out = Struc2Vec(obj,in) %#ok<INUSL>
           out = [];
           if ~isempty(obj.Model), out = [out Struc2Coeff(obj.Model,in)];end
           out = [out in.UserData.Prop];
        end
        function out = Vec2Struc(obj,in) %#ok<INUSL>
                 if obj.nrMC>0
                    modelin = in(1:obj.nrMC);
                    tmpout = Coeff2Struc(obj.Model,modelin);
                    switch class(obj.Model)
                        case 'shapePCA'
                            out.Vertices = tmpout;
                        case 'texturePCA'
                            out.TextureColor = tmpout;
                        case 'textureMapPCA'
                            out.TextureMap = tmpout;
                        case {'appearancePCA' 'ShapeValuePCA'}
                            out = tmpout;
                    end
                    clear tmpout;
                 end
                 if obj.nrP >0
                    out.Prop = in(end-obj.nrP+1:end);                
                 end
        end
        function out = IndexStruc2Vec(obj,in) %#ok<INUSL>
            % To be improved, but probably fine like this
            out = in;
        end
        function out = IndexVec2Struc(obj,in) %#ok<INUSL>
            % To be improved, but probably fine like this
            out = in;
        end
        function out = WeightStruc2Vec(obj,in) %#ok<INUSL>
            % To be improved, but probably fine like this
            out = in;
        end
        function out = WeightVec2Struc(obj,in) %#ok<INUSL>
            % To be improved, but probably fine like this
            out = in;
        end
        function out = getData(obj,in) %#ok<INUSL>
           switch class(in)
                case 'double'
                    out = in;
                case 'struct'
                    if ~isfield(in,'Properties'), error('Properties Field is missing in data structure'); end
                    obj.PropNames = in.Properties.Names;                   
                    out = [];
                    if ~isempty(obj.Model), out = [out obj.Model.Tcoeff'];end
                    size(out)
                    out = [out;in.Properties.Values];
           end
        end
        function [f,scan] = initializeShowPC(obj,coeff)
                 vec = Coeff2Vec(obj,coeff);
                 [f,scan] = initializeShowPC(obj.Model,vec(1:obj.nrMC));
                 disp(num2str(vec(end-obj.nrP+1:end)'));
        end
        function updateShowPC(obj,scan,coeff)
                 vec = Coeff2Vec(obj,coeff);
                 updateShowPC(obj.Model,scan,vec(1:obj.nrMC));
                 disp(num2str(vec(end-obj.nrP+1:end)'));
        end
        function objout = embedInScan(obj,in)
                 checkAverage(obj);
                 objout = clone(obj.Average);
                 struc2obj(objout,in);
                 objout.UserData.Prop = in.Prop;
        end
        function out = getModelCoeff(obj,in)
            out = getVec(obj,in);
            out = out(1:obj.nrMC);
        end
        function out = getShapeCoeff(obj,in)
            switch class(obj.Model)
                case 'shapePCA'
                    out = getModelCoeff(obj,in);
                case {'appearancePCA' 'ShapeValuePCA'}
                    out = getShapeCoeff(obj.Model,getModelCoeff(obj,in));
                otherwise
                    out =[];
            end
        end
        function out = getTextureCoeff(obj,in)
            switch class(obj.Model)
                case {'texturePCA' 'textureMapPCA'}
                    out = getModelCoeff(obj,in);
                case 'appearancePCA'
                    out = getTextureCoeff(obj.Model,getModelCoeff(obj,in));
                case 'ShapeValuePCA'
                    disp('To be implemented');
                otherwise
                    out =[];
            end
        end
        function out = getProperties(obj,in)
            out = getVec(obj,in);
            out = out(obj.nrMC+1:end);
        end
    end
    methods % animate property functions
        function animateProperty(obj,index,path,range,step,coeff,rec)
            % animatePC(obj,pc,range,rec)
            % animate principal component, obj, pc = principal component
            % range [-x x], record true or false;
            if nargin < 7, rec = false; end
            if nargin < 6, coeff = obj.AvgCoeff; end
            if isempty(coeff), coeff = obj.AvgCoeff; end
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
                  MovieFile.fps = 3;
                  MovieFile.compression = 'none';
               end
            end
            coeff = setPropertyValue(obj,coeff,index,range(1),path);
            [f,scan] = initializeShowPC(obj,coeff);
            in = input('Ready? y/n:','s');
            while strcmp(in,'y')
                  MovieFile = runProperty(obj,index,path,range,step,scan,rec,f,MovieFile,coeff);
                  in = input('Again? y/n:','s');
            end
            if rec, MovieFile = close(MovieFile);end %#ok<NASGU>
            if superClass.isH(scan),delete(scan);end
            try delete(f);catch end %#ok<CTCH>
        end
        function movie = runProperty(obj,index,path,range,step,scan,rec,f,movie,coeff)
            %if rec, movie = recordFrame(movie,f);end
            for i=range(1):step:range(2)
                 coeff = setPropertyValue(obj,coeff,index,i,path);
                 updateShowPC(obj,scan,coeff);
                 pause(0.01);
                 if rec,movie = recordFrame(movie,f);end
             end
             for i=1:1:2
                 if rec,movie = recordFrame(movie,f);end
                 pause(0.01);
             end
             for i=range(2):-step:range(1)
                 coeff = setPropertyValue(obj,coeff,index,i,path);
                 updateShowPC(obj,scan,coeff);
                 pause(0.01);
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
    end
end

