classdef BatchMapper < superHandleClass
    % batch processing class, that enables the processing in batch of
    % images to be mapped. The input argument is a the path to a directory,
    % contained pre-defined folders
   properties
       HomeDirectory;% Main directory containing the data and processing steps
       OutputStructure = 'original'% original or list or folder determines the way in which the output files are structured
       Overwrite = false;% OVerwrite existing processed files yes or no
       Preload = true;% Use parallel loop to Preload all scans to indicate and clean (time saving, but memory intensive), this mainly applies to pose indication and scan cleaning;
   end
   properties % Mapping Properties
       TemplateShape;% Shape template to match on all shapes
       TemplatePose;% Pose Points of Template
       TemplateTextureInfo;% Information required to drive the texture mapping from the template
       TemplateReflectIndex;% Give Reflection index if exist to speed up the reflected mapping
       RigidMapper;% RIGID mapper to align the template rigidly with target shapes, after an initial pose aligment
       NonRigidMapper;% NON RIGID mapper tp aling the template non rigiduly with target shapes, after rigid aligment
       GeometricImageRes = 512;
   end
   properties (Dependent = true)
   end
   properties (Hidden = true)
       Status = 'ready';
       MaxN = +inf;
   end
   properties (Dependent = true, Hidden = true)
       pathImages;% path to all IMAGES
       pathRaw;% path to raw data
       path01;% path to IMAGES/01 ORIGINAL IMAGES/
       path11;% path to IMAGES/11 MAT IMAGES/
       path20;% path to IMAGES/20 XML POSE POINTS/
       path21;% path to IMAGES/21 MAT POSE POINTS/
       path30;% path to IMAGES/30 WEM CLEANED/
       path31;% path to IMAGES/31 OBJ CLEANED/
       path32;% path to IMAGES/32 MAT CLEANED/
       path41;% path to IMAGES/41 MAT MAPPED/
       path42;% path to IMAGES/42 TEXTURE MAPS/
       path43;% path to IMAGES/43 MAT REFLECTED MAPPED/
       path51;% path to IMAGES/51 OBJ MAPPED/
       path52;% path to IMAGES/52 OBJ REFLECTED MAPPED/
       path61;% path to IMAGES/61 GEOMETRIC IMAGES/
       path62;% path to IMAGES/62 DEPTH IMAGES/
       pathPar;% path to a temporary folder for parallel computing RAW/PARFOR/
       pathMetaData;% path to metdata
       path00;
   end
   methods % CONSTRUCTOR
       function obj = BatchMapper(varargin)
          obj = obj@superHandleClass(varargin{:});
          obj.RigidMapper = ShapeMapper;
          obj.NonRigidMapper = ShapeMapper;
       end % batch mapper Constructor 
   end
   methods % SETTING
       function obj = set.HomeDirectory(obj,in)
           in = BatchMapper.addBackSlash(in);
           obj.HomeDirectory = in;
           prepHomeDirectory(obj);
       end
       function obj = set.Status(obj,in)
           obj.Status = in;
           disp(obj.Status);
       end
   end
   methods % GETTING
       function out = get.pathImages(obj)
           if isempty(obj.HomeDirectory), out = []; return;end
           out = [obj.HomeDirectory 'IMAGES/'];
       end
       function out = get.pathRaw(obj)
           if isempty(obj.HomeDirectory), out = []; return;end
           out = [obj.HomeDirectory 'RAW/'];
       end
       function out = get.pathPar(obj)
           if isempty(obj.HomeDirectory), out = []; return;end
           out = [obj.HomeDirectory 'RAW/PARFOR/'];
       end
       function out = get.path00(obj)
           if isempty(obj.HomeDirectory), out = []; return;end
           out = [obj.HomeDirectory 'IMAGES/00 BAD IMAGES/'];
       end
       function out = get.path01(obj)
           if isempty(obj.HomeDirectory), out = []; return;end
           out = [obj.HomeDirectory 'IMAGES/01 ORIGINAL IMAGES/'];
       end
       function out = get.path11(obj)
           if isempty(obj.HomeDirectory), out = []; return;end
           out = [obj.HomeDirectory 'IMAGES/11 MAT IMAGES/'];
       end
       function out = get.path20(obj)
           if isempty(obj.HomeDirectory), out = []; return;end
           out = [obj.HomeDirectory 'IMAGES/20 XML POSE POINTS/'];
       end
       function out = get.path21(obj)
           if isempty(obj.HomeDirectory), out = []; return;end
           out = [obj.HomeDirectory 'IMAGES/21 MAT POSE POINTS/'];
       end
       function out = get.path30(obj)
           if isempty(obj.HomeDirectory), out = []; return;end
           out = [obj.HomeDirectory 'IMAGES/30 WEM CLEANED/'];
       end
       function out = get.path31(obj)
           if isempty(obj.HomeDirectory), out = []; return;end
           out = [obj.HomeDirectory 'IMAGES/31 OBJ CLEANED/'];
       end
       function out = get.path32(obj)
           if isempty(obj.HomeDirectory), out = []; return;end
           out = [obj.HomeDirectory 'IMAGES/32 MAT CLEANED/'];
       end
       function out = get.path41(obj)
           if isempty(obj.HomeDirectory), out = []; return;end
           out = [obj.HomeDirectory 'IMAGES/41 MAT MAPPED/'];
       end
       function out = get.path42(obj)
           if isempty(obj.HomeDirectory), out = []; return;end
           out = [obj.HomeDirectory 'IMAGES/42 TEXTURE MAPS/'];
       end
       function out = get.path43(obj)
           if isempty(obj.HomeDirectory), out = []; return;end
           out = [obj.HomeDirectory 'IMAGES/43 MAT REFLECTED MAPPED/'];
       end
       function out = get.path51(obj)
           if isempty(obj.HomeDirectory), out = []; return;end
           out = [obj.HomeDirectory 'IMAGES/51 OBJ MAPPED/'];
       end
       function out = get.path52(obj)
           if isempty(obj.HomeDirectory), out = []; return;end
           out = [obj.HomeDirectory 'IMAGES/52 OBJ REFLECTED MAPPED/'];
       end
       function out = get.path61(obj)
           if isempty(obj.HomeDirectory), out = []; return;end
           out = [obj.HomeDirectory 'IMAGES/61 GEOMETRIC IMAGES/'];
       end
       function out = get.path62(obj)
           if isempty(obj.HomeDirectory), out = []; return;end
           out = [obj.HomeDirectory 'IMAGES/62 DEPTH IMAGES/'];
       end
       function out = get.pathMetaData(obj)
           if isempty(obj.HomeDirectory), out = []; return;end
           out = [obj.HomeDirectory 'METADATA/'];
       end
       function out = get.RigidMapper(obj)
           out = obj.RigidMapper;
           if ~superHandleClass.isH(out), out = []; end
       end
       function out = get.NonRigidMapper(obj)
           out = obj.NonRigidMapper;
           if ~superHandleClass.isH(out), out = []; end
       end
       function out = get.TemplateShape(obj)
           out = obj.TemplateShape;
           if ~superHandleClass.isH(out), out = []; end
       end
       function out = get.TemplatePose(obj)
           out = obj.TemplatePose;
           if ~superHandleClass.isH(out), out = []; end
       end    
   end
   methods % MAPPING STEPS
       function [info,inputfiles] = step1Convert2MAT(obj)
           inputpath = obj.path01;outputpath = obj.path11;
           [inputfiles,outputfiles,done] = prepIO(obj,inputpath,outputpath,'.obj','.mat');
           if isempty(inputfiles), info = []; return; end% nothing to process
           info = ones(1,length(done));index = find(done==0);N = length(index);
           if N==0, return; end% nothing to process
           shortinputfiles = inputfiles(index);shortoutputfiles = outputfiles(index);
           shortinfo = ones(1,N);
           shape = shape3D;
           [parforpath,ID] = openParfor(obj,N);
           parfor i=1:N
              try
                  forshape = clone(shape);
                  importWavefront(forshape,shortinputfiles{i}.FileName,shortinputfiles{i}.FullFolder,[]);
                  rgb = forshape.VertexRGB;forshape.TextureMap = [];forshape.VertexRGB = rgb;
                  cleanManifold(shape);% very frequently nan vertices or triangles referencing to non-vertices exist and are removed here
                  % The idea is not to drag the texture maps along, as these are memory intensive
                  % Hence when doing texturemapping, the textures are to be
                  % gathered from the original folder again. Colors per vertex
                  % are created and used during the subsequent steps including
                  % landmarking
                  BatchMapper.saveMatFile(shortoutputfiles{i},forshape);
              catch
                  shortinfo(i) = 0;
              end
              BatchMapper.progressParfor; 
           end
           closeParfor(obj,parforpath,ID);
           info(index) = shortinfo;
       end
       function [info,inputfiles] = step2PoseIndication(obj,inputpath)
           if nargin<2, inputpath = obj.path11;end% default input path is from converted images, but can be set differently
           outputpath = obj.path21;% define input <-> output folders
           [inputfiles,outputfiles,done] = prepIO(obj,inputpath,outputpath,'.mat','.mat');% get input and output files
           if isempty(inputfiles), info = []; return; end% nothing to process
           info = ones(1,length(done));index = find(done==0);N = length(index);
           if N==0, return; end% nothing to process
           shortinputfiles = inputfiles(index);shortoutputfiles = outputfiles(index);shortinfo = ones(1,N);
           [shapes] = preLoader(obj,shortinputfiles);% preload shapes
           [parforpath,ID] = openParfor(obj,N);
           for i=1:1:N
               if ~isempty(shapes)
                  shape = shapes{i};
               else
                  shape = BatchMapper.loadMatFile(shortinputfiles{i});
               end
               if~isempty(shape.VertexRGB), shape.ColorMode = 'texture';else, shape.ColorMode = 'single';end
               shape.ViewMode = 'solid';
               v = viewer(shape);
               v.SelectionMode = 'landmark';
               v.SelectionActive = true;
               v.Figure.Position = [98.1429   14.1250  155.5714   51.0625];
               v.Tag = shortinputfiles{i}.Name;
               if strcmp(shape.ColorMode,'single'), v.SceneLightVisible = true;v.SceneLightLinked = true;end
               waitfor(v);
               % perform a check(s)
               if isfield(shape.UserData,'LandmarkSelection'),posepoints = clone(shape.UserData.LandmarkSelection);else, shortinfo(i) = 0; BatchMapper.progressParfor; continue; end
               if posepoints.nVertices==5
                  BatchMapper.saveMatFile(shortoutputfiles{i},posepoints);
               else
                   disp('Warning, failed to indicate 5 pose landmarks, results are not saved');
                   shortinfo(i) = 0;
               end
              BatchMapper.progressParfor;
           end
           closeParfor(obj,parforpath,ID);
           info(index) = shortinfo;
       end
       function [info,inputfiles] = step3Cleaning(obj,inputpath)
           if nargin<2, inputpath = obj.path11;end% default input path is from converted images, but can be set differently
           outputpath = obj.path32;% define input <-> output folders
           [inputfiles,outputfiles,done] = prepIO(obj,inputpath,outputpath,'.mat','.mat');% get input and output files
           if isempty(inputfiles), info = []; return; end% nothing to process
           info = ones(1,length(done));index = find(done==0);N = length(index);
           if N==0, return; end% nothing to process
           shortinputfiles = inputfiles(index);shortoutputfiles = outputfiles(index);shortinfo = ones(1,N);
           [shapes] = preLoader(obj,shortinputfiles);% preload shapes
           [parforpath,ID] = openParfor(obj,N);
           for i=1:1:N
               if ~isempty(shapes)
                  shape = shapes{i};
               else
                  shape = BatchMapper.loadMatFile(shortinputfiles{i});
               end
               if~isempty(shape.VertexRGB), shape.ColorMode = 'texture';end
               shape.ViewMode = 'wireframe';
               v = viewer(shape);
               v.SelectionMode = 'brush';
               v.SelectionActive = true;
               v.Figure.Position = [98.1429   14.1250  155.5714   51.0625];
               v.Tag = shortinputfiles{i}.Name;
               waitfor(v);
               cleanManifold(shape);% removing nan and non-connected vertices
               BatchMapper.saveMatFile(shortoutputfiles{i},shape);
               BatchMapper.progressParfor;
           end
           closeParfor(obj,parforpath,ID);
           info(index) = shortinfo;
       end
       function [info,inputfiles] = step4MapShape(obj,reflect)
           if nargin<2, reflect = false; end% reflect on creates output for relfected original shapes
           inputpath = obj.path32;posepath = obj.path21;
           if ~reflect, outputpath = obj.path41; else, outputpath = obj.path43;end
           [inputfiles,outputfiles,done] = prepIO(obj,inputpath,outputpath,'.mat','.mat');% get input and output files
           if isempty(inputfiles), info = []; return; end% nothing to process
           info = ones(1,length(done));index = find(done==0);N = length(index);
           if N==0, return; end% nothing to process
           shortinputfiles = inputfiles(index);shortoutputfiles = outputfiles(index);shortinfo = ones(1,N);
           % cloning template information this not to push the complete object into the parfor loop
           template = clone(obj.TemplateShape);template.UserData = [];template.TextureMap = [];% traveling light in the parfor loop
           templatepose = clone(obj.TemplatePose);RM = obj.RigidMapper;NRM = obj.NonRigidMapper;
           normalOrientationTemplate = getNormalOrientation(template);
           [parforpath,ID] = openParfor(obj,N);
           parfor i=1:N
              shape = BatchMapper.loadMatFile(shortinputfiles{i});cleanManifold(shape);% load and clean
              if ~getNormalOrientation(shape)==normalOrientationTemplate, shape.FlipNormals = true; end% check the orientation of the normals matches
              posefile = clone(shortinputfiles{i});posefile.MainFolder = posepath; % Retrieving Pose Points to start the mapping
              if ~fileExist(posefile)% try to see if you can find the posefile directly or in another subfolder
                 folder = BatchMapper.locateFile(posepath,posefile.FileName);% see if the file does not exist in another subfolder
                 if isempty(folder), shortinfo(i) = 2,BatchMapper.progressParfor;continue;end% no existing pose points, mapping cannot be initiated
                 posefile.MainFolder = folder;posefile.SubFolder = [];% correct the location of the file (which appeared to be in a different subfolder than the inputfile
              end
              shapepose = BatchMapper.loadMatFile(posefile);% this is on the expected location
              try
                 if reflect, [shape,shapepose] = BatchMapper.reflectShape(shape,shapepose);end 
                 mappedshape = BatchMapper.mapShape(shape,shapepose,template,templatepose,RM,NRM);
                 BatchMapper.saveMatFile(shortoutputfiles{i},mappedshape);
              catch
                  shortinfo(i) = 0;% fialed to map, requires further investigation
              end
              BatchMapper.progressParfor; 
           end
           closeParfor(obj,parforpath,ID);
           info(index) = shortinfo;
       end
       function [info,inputfiles] = step4MapTexture(obj)
           inputpath = obj.path41;outputpath = obj.path42;texturepath = obj.path01;
           %origpath = obj.path32;
           origpath = obj.path11;
           [inputfiles,outputfiles,done] = prepIO(obj,inputpath,outputpath,'.mat','.bmp');% get input and output files
           if isempty(inputfiles), info = []; return; end% nothing to process
           info = ones(1,length(done));index = find(done==0);N = length(index);
           if N==0, return; end% nothing to process
           shortinputfiles = inputfiles(index);shortoutputfiles = outputfiles(index);shortinfo = ones(1,N);
           % cloning template information this not to push the complete object into the parfor loop
           template = clone(obj.TemplateShape);% traveling light in the parfor loop
           textureinfo = obj.TemplateTextureInfo;
           [parforpath,ID] = openParfor(obj,N);
           parfor i=1:N
              shape = BatchMapper.loadMatFile(shortinputfiles{i});cleanManifold(shape);% load mapped shape result
              origfile = clone(shortinputfiles{i});origfile.MainFolder = origpath; % Retrieving Original Shape
              if ~fileExist(origfile)% try to see if you can find the original file directly or in another subfolder
                 folder = BatchMapper.locateFile(origpath,origfile.FileName);% see if the file does not exist in another subfolder
                 if isempty(folder), shortinfo(i) = 2,BatchMapper.progressParfor;continue;end% no existing original file, texture mapping cannot be initiated
                 origfile.MainFolder = folder;origfile.SubFolder = [];% correct the location of the file (which appeared to be in a different subfolder than the inputfile
              end
              origshape = BatchMapper.loadMatFile(origfile);% this is on the expected location
              texturefile = clone(shortinputfiles{i});texturefile.MainFolder = texturepath; texturefile.Extension = '.bmp';% Retrieving Original TextureMap
              if ~fileExist(texturefile)% try to see if you can find the texture file directly or in another subfolder
                 folder = BatchMapper.locateFile(texturepath,texturefile.FileName);% see if the file does not exist in another subfolder
                 if isempty(folder), texturefile.Extension = '.png'; folder = BatchMapper.locateFile(texturepath,texturefile.FileName); end % change extension to .png
                 if isempty(folder), shortinfo(i) = 3,BatchMapper.progressParfor;continue;end% no existing texture file, texture mapping cannot be initiated
                 texturefile.MainFolder = folder;texturefile.SubFolder = [];% correct the location of the file (which appeared to be in a different subfolder than the inputfile
              end
              try
                origshape.TextureMap = imread([texturefile.FullName]);% it is sometimes possible that the texture file is corrupted!  
                mapTexture(shape,template,textureinfo,origshape);
                BatchMapper.saveMatFile(shortinputfiles{i},shape);
                imwrite(shape.TextureMap,shortoutputfiles{i}.FullName,'bmp');
              catch
                shortinfo(i) = 0;% fialed to map, requires further investigation  
              end
              BatchMapper.progressParfor; 
           end
           closeParfor(obj,parforpath,ID);
           info(index) = shortinfo;
       end
       function [info,inputfiles] = step5ExportWavefront(obj,reflect)
           if nargin<2, reflect = false; end% reflect on creates output for relfected original shapes
           texturepath = obj.path42;
           if ~reflect, inputpath = obj.path41; else, inputpath = obj.path43;end
           if ~reflect, outputpath = obj.path51; else, outputpath = obj.path52;end
           [inputfiles,outputfiles,done] = prepIO(obj,inputpath,outputpath,'.mat','.obj');% get input and output files
           if isempty(inputfiles), info = []; return; end% nothing to process
           info = ones(1,length(done));index = find(done==0);N = length(index);
           if N==0, return; end% nothing to process
           shortinputfiles = inputfiles(index);shortoutputfiles = outputfiles(index);shortinfo = ones(1,N);
           [parforpath,ID] = openParfor(obj,N);
           parfor i=1:N
              try
                  shape = BatchMapper.loadMatFile(shortinputfiles{i});cleanManifold(shape);% load and clean
                  texturefile = clone(shortinputfiles{i});texturefile.MainFolder = texturepath; texturefile.Extension = '.bmp';% Retrieving Original TextureMap
                  if ~fileExist(texturefile)% try to see if you can find the texture file directly or in another subfolder
                     folder = BatchMapper.locateFile(texturepath,texturefile.FileName);% see if the file does not exist in another subfolder
                     texturefile.MainFolder = folder;texturefile.SubFolder = [];% correct the location of the file (which appeared to be in a different subfolder than the inputfile
                  end
                  if ~isempty(texturefile.MainFolder),shape.TextureMap = imread([texturefile.FullName]);end
                  exportWavefront(shape,shortoutputfiles{i}.Name,shortoutputfiles{i}.FullFolder);
              catch
                   shortinfo(i) = 0;% fialed to map, requires further investigation
              end
              BatchMapper.progressParfor; 
           end
           closeParfor(obj,parforpath,ID);
           info(index) = shortinfo;
       end
       function [info,inputfiles] = step6GeometricImage(obj)
           inputpath = obj.path41;outputpath = obj.path61;
           [inputfiles,outputfiles,done] = prepIO(obj,inputpath,outputpath,'.mat','.png');% get input and output files
           if isempty(inputfiles), info = []; return; end% nothing to process
           info = ones(1,length(done));index = find(done==0);N = length(index);
           if N==0, return; end% nothing to process
           shortinputfiles = inputfiles(index);shortoutputfiles = outputfiles(index);shortinfo = ones(1,N);
           % cloning template information this not to push the complete object into the parfor loop
           template = clone(obj.TemplateShape);
           range(:,1) = min(template.Vertices);
           range(:,2) = max(template.Vertices);
           range = range.*1.2;
           res = obj.GeometricImageRes;
           [parforpath,ID] = openParfor(obj,N);
           parfor i=1:N
               try
                  shape = BatchMapper.loadMatFile(shortinputfiles{i});cleanManifold(shape);% load and clean
                  shape = BatchMapper.alignWithTemplate(shape,template);
                  out = geometricImage(shape,res,range);
                  imwrite(out,shortoutputfiles{i}.FullName,'png');
               catch
                  shortinfo(i) = 0;% fialed to map, requires further investigation 
               end
              BatchMapper.progressParfor; 
           end
           closeParfor(obj,parforpath,ID);
           info(index) = shortinfo;
       end
       function [info,inputfiles] = step6DepthImage(obj)
           inputpath = obj.path41;outputpath = obj.path62;
           [inputfiles,outputfiles,done] = prepIO(obj,inputpath,outputpath,'.mat','.png');% get input and output files
           if isempty(inputfiles), info = []; return; end% nothing to process
           info = ones(1,length(done));index = find(done==0);N = length(index);
           if N==0, return; end% nothing to process
           shortinputfiles = inputfiles(index);shortoutputfiles = outputfiles(index);shortinfo = ones(1,N);
           % cloning template information this not to push the complete object into the parfor loop
           template = clone(obj.TemplateShape);
           range(:,1) = min(template.Vertices);
           range(:,2) = max(template.Vertices);
           range = range.*1.2;
           res = obj.GeometricImageRes;
           [parforpath,ID] = openParfor(obj,N);
           parfor i=1:N
               try
                  shape = BatchMapper.loadMatFile(shortinputfiles{i});cleanManifold(shape);% load and clean
                  shape = BatchMapper.alignWithTemplate(shape,template);
                  out = geometricImage(shape,res,range);
                  out = squeeze(out(:,:,3));
                  imwrite(out,shortoutputfiles{i}.FullName,'png');
               catch
                  shortinfo(i) = 0;% fialed to map, requires further investigation 
               end
              BatchMapper.progressParfor; 
           end
           closeParfor(obj,parforpath,ID);
           info(index) = shortinfo;
       end
   end
   methods % INTERFACING
       function prepHomeDirectory(obj)
           BatchMapper.prepFolder(obj.HomeDirectory);
           BatchMapper.prepFolder(obj.pathImages);
           BatchMapper.prepFolder(obj.pathRaw);
           BatchMapper.prepFolder(obj.pathPar);
           BatchMapper.prepFolder(obj.pathMetaData);
           BatchMapper.prepFolder(obj.path00);
           BatchMapper.prepFolder(obj.path01);
           BatchMapper.prepFolder(obj.path11);
           %BatchMapper.prepFolder(obj.path20);
           BatchMapper.prepFolder(obj.path21);
           %BatchMapper.prepFolder(obj.path30);
           %BatchMapper.prepFolder(obj.path31);
           BatchMapper.prepFolder(obj.path32);
           BatchMapper.prepFolder(obj.path41);
           BatchMapper.prepFolder(obj.path42);
           BatchMapper.prepFolder(obj.path43);
           BatchMapper.prepFolder(obj.path51);
           BatchMapper.prepFolder(obj.path52);
           BatchMapper.prepFolder(obj.path61);
           BatchMapper.prepFolder(obj.path62);
           cd(obj.HomeDirectory);
       end
       function [inputfiles,outputfiles,done] = prepIO(obj,inputpath,outputpath,inputformat,outputformat)
                disp(' ');disp('...   preparing IO   ...');disp(' ');
                % Initialization
                 outputpath = BatchMapper.addBackSlash(outputpath);
                 inputpath = BatchMapper.addBackSlash(inputpath);
                % STEP 1 ORGANIZING INPUT   
                 [inputfiles,ninput] = BatchMapper.getFiles(inputpath,inputformat);
                 if ninput==0, disp('nothing to process'); outputfiles = []; done = []; return; end
                % STEP2 CHECK ON AND ORGANIZE EXISTING FILES
                 [existingfiles,nexisting] = BatchMapper.getFiles(outputpath,outputformat);
                 if nexisting>0
                    existingarray = cell(1,nexisting);
                    for i=1:nexisting
                        existingarray{i}=existingfiles{i}.Name;
                    end
                 else
                    existingarray = {};
                 end 
                % STEP2 ORGANIZING OUTPUT
                ninput = min(ninput,obj.MaxN);% you can restrict yourself in numbers by setting MaxN to a maximum number of files to process, this is usefull when developping subroutines
                inputfiles = inputfiles(1:ninput);
                outputfiles = cell(1,ninput);
                structure = lower(obj.OutputStructure);
                done = zeros(1,ninput);
                [path,parID] = openParfor(obj,ninput);
                parfor i=1:1:ninput
                    %disp(num2str(i));
                    % create designated outputfile name and folder
                    outputfiles{i} = clone(inputfiles{i});
                    outputfiles{i}.Extension = outputformat;
                    outputfiles{i}.MainFolder = outputpath;
                    switch structure  
                        case 'original'
                            % nothing to do, cloned the inputfiles's subfolder
                        case 'list'
                            outputfiles{i}.SubFolder = [];% all files into a single output folder
                        case 'folder'
                            outputfiles{i}.SubFolder = [outputfiles{i}.Name '/'];% all files, each in his own folder  
                    end
                    BatchMapper.prepFolder(outputfiles{i}.FullFolder);% prepare the designated outputfolder
                    index = find(strcmp(outputfiles{i}.Name,existingarray));% see if an existing file exists
                    if isempty(index), BatchMapper.progressParfor; continue;end% if not then simply continue;
                    if strcmp(existingfiles{index(1)}.FullFolder,outputfiles{i}.FullFolder), done(i) = true; BatchMapper.progressParfor; continue; end%#ok<*PFBNS> % file exists, code true, and is in de right folder
                    % file exist, code = 1, but is in the wrong folder, so needs to be moved.
                    succes = movefile(existingfiles{index(1)}.FullName,outputfiles{i}.FullName,'f');
                    if succes, done(i) = true; end
                    BatchMapper.progressParfor;
                end
                if obj.Overwrite, done = zeros(1,ninput); end
                closeParfor(obj,path,parID)
                disp(' ');
                disp(['...   Ready to process: ' num2str(ninput-sum(done)) ' of ' num2str(ninput) ' Files   ...']);% how many to process?
                disp('  ');
       end
       function [path,ID] = openParfor(obj,N)
           path = pwd;
           cd(obj.pathPar);
           ID = num2str(randi(10000,1));
           BatchMapper.prepFolder(ID);cd(ID);
           %mkdir(ID);cd(ID);
           parfor_progress(N);
       end
       function closeParfor(obj,path,ID)
           cd([obj.pathPar ID]);
           parfor_progress(0);
           cd(path);
           try% this fails under windows for some reason
            rmdir([obj.pathPar ID '/'],'s');
           catch
           end
       end
       function presetMapperSettings(obj,type)
           if isempty(obj.RigidMapper), obj.RigidMapper = ShapeMapper; end
           if isempty(obj.NonRigidMapper), obj.NonRigidMapper = ShapeMapper; end
           RM = obj.RigidMapper;NRM = obj.NonRigidMapper;% shorter notation througout function
           RM.Display = false;NRM.Display = false;
           switch lower(type)
               case 'faces'
                   % RIGID MAPPER
                   RM.NumIterations = 30;
                   RM.InlierKappa = 3;
                   RM.TransformationType = 'rigid';
                   RM.UseScaling = true;
                   RM.CorrespondencesNumNeighbours = 3;
                   RM.CorrespondencesFlagThreshold = 0.9;
                   RM.CorrespondencesSymmetric = true;
                   RM.CorrespondencesEqualizePushPull = false;
                   RM.InlierUseOrientation = true;
                   
                   RM.FlagFloatingBoundary = true;
                   RM.FlagTargetBoundary = true;
                   
                   RM.FlagTargetBadlySizedTriangles = true;
                   RM.TriangleSizeZscore = 6;
                   
                   RM.UpSampleTarget = false;
                   
                   % NON RIGID MAPPER
                   nr = 200;
                   NRM.NumIterations = nr;
                   NRM.CorrespondencesSymmetric = true;
                   NRM.CorrespondencesNumNeighbours = 3;
                   NRM.CorrespondencesFlagThreshold = 0.9;
                   NRM.CorrespondencesUseOrientation = true;
                   NRM.CorrespondencesEqualizePushPull =false;
                   NRM.InlierKappa = 12;
                   NRM.InlierUseOrientation = true;
                   NRM.FlagFloatingBoundary = true;
                   NRM.FlagTargetBoundary = true;
                   
                   NRM.FlagTargetBadlySizedTriangles = true;
                   NRM.TriangleSizeZscore = 6;
                   
                   NRM.UpSampleTarget = false;
                   
                   NRM.UseScaling = 1;
                   NRM.TransformationType = 'nonrigid';
                   NRM.TransformSigma = 3;
                   NRM.TransformNumViscousIterationsStart = nr;
                   NRM.TransformNumViscousIterationsEnd = 1;
                   NRM.TransformNumElasticIterationsStart = nr;
                   NRM.TransformNumElasticIterationsEnd = 1;
                   NRM.TransformNumNeighbors = 80;
               case 'irregular meshes'
                   % RIGID MAPPER
                   RM.NumIterations = 30;
                   RM.InlierKappa = 3;
                   RM.TransformationType = 'rigid';
                   RM.UseScaling = true;
                   RM.CorrespondencesNumNeighbours = 3;
                   RM.CorrespondencesFlagThreshold = 0.9;
                   RM.CorrespondencesSymmetric = true;
                   RM.CorrespondencesEqualizePushPull = false;
                   RM.InlierUseOrientation = true;
                   
                   RM.FlagFloatingBoundary = true;
                   RM.FlagTargetBoundary = true;
                   
                   RM.FlagTargetBadlySizedTriangles = false;
                   RM.TriangleSizeZscore = 6;
                   
                   RM.UpSampleTarget = true;
                   RM.UpSampleMode = 'size';
                   RM.UpSampleVal = 1.5;
                   % NON RIGID MAPPER
                   nr = 200;
                   NRM.NumIterations = nr;
                   NRM.CorrespondencesSymmetric = true;
                   NRM.CorrespondencesNumNeighbours = 3;
                   NRM.CorrespondencesFlagThreshold = 0.9;
                   NRM.CorrespondencesUseOrientation = true;
                   NRM.CorrespondencesEqualizePushPull =false;
                   NRM.InlierKappa = 12;
                   NRM.InlierUseOrientation = true;
                   
                   NRM.FlagFloatingBoundary = true;
                   NRM.FlagTargetBoundary = true;
                   
                   NRM.FlagTargetBadlySizedTriangles = false;
                   NRM.TriangleSizeZscore = 6;
                   % this is the key change
                   NRM.UpSampleTarget = true;
                   NRM.UpSampleMode = 'size';
                   NRM.UpSampleVal = 1.5;
                   
                   NRM.UseScaling = 1;
                   NRM.TransformationType = 'nonrigid';
                   NRM.TransformSigma = 3;
                   NRM.TransformNumViscousIterationsStart = nr;
                   NRM.TransformNumViscousIterationsEnd = 1;
                   NRM.TransformNumElasticIterationsStart = nr;
                   NRM.TransformNumElasticIterationsEnd = 1;
                   NRM.TransformNumNeighbors = 80;    
               otherwise
           end
           
       end
       function [shapes] = preLoader(obj,inputfiles)
           N = length(inputfiles);
           if obj.Preload, shapes = cell(1,N); else, shapes = []; return; end
           % Do a parallel preload to gain speed
           disp(' ');disp('...   Loading   ...');disp(' ');
           [parforpath,ID] = openParfor(obj,N);
           parfor i=1:N
              shapes{i} = BatchMapper.loadMatFile(inputfiles{i});
              BatchMapper.progressParfor;
           end
           closeParfor(obj,parforpath,ID);
           disp(' ');disp('...   Loaded   ...');disp(' ');
       end
       function out = getStats(obj)
           [~,out(1)] = BatchMapper.getFiles(obj.path01,'.obj');
           disp(['# original :' num2str(out(1))]);
           [~,out(2)] = BatchMapper.getFiles(obj.path11,'.mat');
           disp(['# converted :' num2str(out(2))]);
           [~,out(3)] = BatchMapper.getFiles(obj.path21,'.mat');
           disp(['# indicated :' num2str(out(3))]);
           [~,out(4)] = BatchMapper.getFiles(obj.path32,'.mat');
           disp(['# cleaned:' num2str(out(4))]);
           [~,out(5)] = BatchMapper.getFiles(obj.path41,'.mat');
           disp(['# mapped :' num2str(out(5))]);
           [~,out(6)] = BatchMapper.getFiles(obj.path42,'.bmp');
           disp(['# texture :' num2str(out(6))]);
       end
       function [out,cellout] = getProcessingReport(obj,fileid)
           %cellout = [];
           headers = {'ID' 'Original' 'Converted' 'Pose' 'Cleaned' 'Mapped' 'Texture' 'Geometric2D' 'Depth2D'};
           if nargin<2
              [inputfiles,nfiles] = BatchMapper.getFiles(obj.path01,'.obj');
              fileid = cell(1,nfiles);
              for i=1:1:nfiles
                  fileid{i} = inputfiles{i}.Name;
              end
           end
           if iscell(fileid)
              nfiles = length(fileid);
              out = ones(nfiles,8);
              [parforpath,ID] = openParfor(obj,nfiles);
              %parfor i=1:nfiles
              for i=1:nfiles
                 out(i,:) = getProcessingReport(obj,fileid{i});
                 BatchMapper.progressParfor;
              end
              closeParfor(obj,parforpath,ID);
              cellout = [fileid(:) num2cell(out)];
              cellout = [headers;cellout];
              return;
           end
           out = ones(1,8);
           % easy implementation, but might be slow for very large datasets
           %original file
           folder = BatchMapper.locateFile(obj.path01,[fileid '.obj']);
           if isempty(folder), out(1) = 0;end
           %converted file
           folder = BatchMapper.locateFile(obj.path11,[fileid '.mat']);
           if isempty(folder), out(2) = 0;end
           % pose landmarks in mat file
           folder = BatchMapper.locateFile(obj.path21,[fileid '.mat']);
           if isempty(folder), out(3) = 0;end
           % cleaned mat file
           folder = BatchMapper.locateFile(obj.path32,[fileid '.mat']);
           if isempty(folder), out(4) = 0;end
           % mapped mat file
           folder = BatchMapper.locateFile(obj.path41,[fileid '.mat']);
           if isempty(folder), out(5) = 0;end
           % texture mapping
           folder = BatchMapper.locateFile(obj.path42,[fileid '.bmp']);
           if isempty(folder), out(6) = 0;end
           % geometric image
           folder = BatchMapper.locateFile(obj.path61,[fileid '.png']);
           if isempty(folder), out(7) = 0;end
           % depth image
           folder = BatchMapper.locateFile(obj.path62,[fileid '.png']);
           if isempty(folder), out(8) = 0;end
           cellout = [{fileid} num2cell(out)];
           cellout = [headers;cellout];
%            % ORIGINAL FILE
%            file = BatchFile;
%            file.FileName = [fileid '.obj'];
%            file.MainFolder = obj.path01;
%            folder = BatchMapper.locateFile(file.MainFolder,file.FileName);
%            if isempty(folder)
%               out(1) = 0;
%            end
%            if length(folder)>length(file.MainFolder)
%               file.SubFolder = folder(length(file.MainFolder)+1:end); 
%            end
%            folder = BatchMapper.addBackSlash(folder);
%            delete([folder fileid ext]);
           
           
           
%            [inputfiles,nfiles] = BatchMapper.getFiles(obj.path01,'.obj');
%            if nargin<2,fileid end
%            Names = cell(1,nfiles);
%            
           
           
       end
       function removeFile(obj,fileid,path)
           switch lower(path)
               case 'converted'
                   path = obj.path11;
                   ext = '.mat';
               case 'pose'
                   path = obj.path21;
                   ext = '.mat';
               case 'cleaned'
                   path = obj.path32;
                   ext = '.mat';
               case 'mapped'
                   path = obj.path41;
                   ext = '.mat';
               case 'texture'
                   path = obj.path42;
                   ext = '.bmp';
               case 'geometric'
                   path = obj.path61;
                   ext = '.png';
               case 'depth'
                   path = obj.path62;
                   ext = '.png';    
               case 'all'
                   removeFile(obj,fileid,'converted');
                   removeFile(obj,fileid,'pose');
                   removeFile(obj,fileid,'cleaned');
                   removeFile(obj,fileid,'processed');
                   return;
               case 'processed'
                   removeFile(obj,fileid,'mapped');
                   removeFile(obj,fileid,'texture');
                   removeFile(obj,fileid,'geometric');
                   removeFile(obj,fileid,'depth');
                   return;
           end
           folder = BatchMapper.locateFile(path,[fileid ext]);
           if isempty(folder), return; end
           folder = BatchMapper.addBackSlash(folder);
           delete([folder fileid ext]);
       end
       function [v,cleanedshape,mappedshape] = inspectMappedResult(obj,fileid)
           mappedpath = obj.path41;
           cleanedpath = obj.path32;
           inputfile = BatchFile;
           inputfile.FileName = [fileid '.mat'];
           inputfile.MainFolder = mappedpath;
           folder = BatchMapper.locateFile(inputfile.MainFolder,inputfile.FileName);
           if isempty(folder), disp('File ID not found in mapped results');return;end
           if length(folder)>length(inputfile.MainFolder)
              inputfile.SubFolder = folder(length(inputfile.MainFolder)+1:end);
           end
           mappedshape = BatchMapper.loadMatFile(inputfile);
           if ~isempty(mappedshape.TextureMap), figure;imshow(mappedshape.TextureMap); end
           
           v = viewer(mappedshape);
           if isfield(mappedshape.UserData,'InlierWeights')
              mappedshape.VertexValue = mappedshape.UserData.InlierWeights;
              mappedshape.ColorMode = 'Indexed';
           end
           cleanedfile = clone(inputfile);cleanedfile.MainFolder = cleanedpath; % Retrieving Original cleaned Shape
           if ~fileExist(cleanedfile)% try to see if you can find the original file directly or in another subfolder
              folder = BatchMapper.locateFile(cleanedpath,cleanedfile.FileName);% see if the file does not exist in another subfolder
              if isempty(folder), disp('Cleaned file could not be found'); return;end% no existing original file, texture mapping cannot be initiated
              cleanedfile.MainFolder = folder;cleanedfile.SubFolder = [];% correct the location of the file (which appeared to be in a different subfolder than the inputfile
           end
           cleanedshape = BatchMapper.loadMatFile(cleanedfile);
           viewer(cleanedshape,v);
           
       end
       function [mappedshapenew,mappedshape] = illustrateMapping(obj,fileid)
           mappedpath = obj.path41;
           cleanedpath = obj.path32;
           posepath = obj.path21;
           
           inputfile = BatchFile;
           inputfile.FileName = [fileid '.mat'];
           inputfile.MainFolder = mappedpath;
           folder = BatchMapper.locateFile(inputfile.MainFolder,inputfile.FileName);
           if isempty(folder), disp('File ID not found in mapped results');return;end
           if length(folder)>length(inputfile.MainFolder)
              inputfile.SubFolder = folder(length(inputfile.MainFolder)+1:end);
           end
           try
            mappedshape = BatchMapper.loadMatFile(inputfile);
           catch
            mappedshape = [];
           end
           
           cleanedfile = clone(inputfile);cleanedfile.MainFolder = cleanedpath; % Retrieving Original cleaned Shape
           if ~fileExist(cleanedfile)% try to see if you can find the original file directly or in another subfolder
              folder = BatchMapper.locateFile(cleanedpath,cleanedfile.FileName);% see if the file does not exist in another subfolder
              if isempty(folder), disp('Cleaned file could not be found'); return;end% no existing original file, texture mapping cannot be initiated
              cleanedfile.MainFolder = folder;cleanedfile.SubFolder = [];% correct the location of the file (which appeared to be in a different subfolder than the inputfile
           end
           cleanedshape = BatchMapper.loadMatFile(cleanedfile);
           
           posefile = clone(inputfile);posefile.MainFolder = posepath; % Retrieving Original cleaned Shape
           if ~fileExist(posefile)% try to see if you can find the original file directly or in another subfolder
              folder = BatchMapper.locateFile(posepath,posefile.FileName);% see if the file does not exist in another subfolder
              if isempty(folder), disp('Cleaned file could not be found'); return;end% no existing original file, texture mapping cannot be initiated
              posefile.MainFolder = folder;posefile.SubFolder = [];% correct the location of the file (which appeared to be in a different subfolder than the inputfile
           end
           poseshape = BatchMapper.loadMatFile(posefile);
           
           template = clone(obj.TemplateShape);template.UserData = [];template.TextureMap = [];% traveling light in the parfor loop
           templatepose = clone(obj.TemplatePose);RM = obj.RigidMapper;NRM = obj.NonRigidMapper;
           normalOrientationTemplate = getNormalOrientation(template);
           
%            templateresolution = median(edgeLengths(template));
%            shaperesolution = median(edgeLengths(cleanedshape));
%            
% %            NRM.NumIterations = 100;
% %            RM.InlierUseOrientation = false;
% %            RM.InlierKappa = 12;
           if ~getNormalOrientation(cleanedshape)==normalOrientationTemplate, cleanedshape.FlipNormals = true; end
           mappedshapenew = BatchMapper.mapShape(cleanedshape,poseshape,template,templatepose,RM,NRM,true);
           
       end
       function chmodImages(obj)
           cd(obj.HomeDirectory);
           if ispc, return;end% you are dealing with a windows computer, no permissions to change
           unix('chmod -R a+rwx IMAGES 2> /dev/null');
           unix('chmod -R a+rwx RAW 2> /dev/null');
       end
       function removeBadImages(obj,fileid)
           if iscell(fileid)
              nfiles = length(fileid);
              for i=1:nfiles
                 disp(['removing image ' num2str(i) ' from ' num2str(nfiles) ' images...']);
                 removeBadImages(obj,fileid{i}); 
              end
              return;
           end
           % locate and move original files
           folder = BatchMapper.locateFile(obj.path01,[fileid '.obj']);
           if ~isempty(folder)
              try movefile([folder fileid '.obj'],[obj.path00 fileid '.obj']);catch, end
              try movefile([folder fileid '.mtl'],[obj.path00 fileid '.mtl']);catch, end
              try movefile([folder fileid '.bmp'],[obj.path00 fileid '.bmp']);catch, end
              try movefile([folder fileid '.png'],[obj.path00 fileid '.png']);catch, end
           end
           % remove any converted, indicated, cleaned and processed files;
           removeFile(obj,fileid,'all');
       end
       function delete(obj)
            chmodImages(obj);
       end
   end
   methods % OUTPUT INTERFACING
       function [Names,Shapes,PercOutlier,ReflectedShapes] = loadMappedData(obj)
           inputpath = obj.path41;% path to the mapped images;
           [inputfiles,nfiles] = BatchMapper.getFiles(inputpath,'.mat');
           Names = cell(1,nfiles);
           Shapes = nan*zeros(obj.TemplateShape.nVertices,3,nfiles);
           if nargout==4, ReflectedShapes = Shapes; reflect = true; else, ReflectedShapes = []; reflect = false; end
           [parforpath,ID] = openParfor(obj,nfiles);
           reflectindex = obj.TemplateReflectIndex;
           PercOutlier = nan*zeros(1,nfiles);
           parfor i=1:nfiles
               shape = BatchMapper.loadMatFile(inputfiles{i});
               if isfield(shape.UserData,'InlierWeights')
                  val = shape.UserData.InlierWeights;
                  PercOutlier(i) = (sum(val)/length(val))*100; 
               end
               Names{i} = inputfiles{i}.Name;
               Shapes(:,:,i) = shape.Vertices;
               if reflect
                   newshape = reflectShape(shape,reflectindex);
                   ReflectedShapes(:,:,i) = newshape.Vertices;
               end
               BatchMapper.progressParfor; 
           end
           closeParfor(obj,parforpath,ID);
       end
   end
   methods % CONVERTORS (mainly from old MevisLab data steps, and meshObj from Penn State)
       function [info,inputfiles] = Xml2MatPoseFiles(obj,inputpath)
           if nargin<2, inputpath = obj.path20;end
           outputpath = obj.path21;
           [inputfiles,outputfiles,done] = prepIO(obj,inputpath,outputpath,'.xml','.mat');
           if isempty(inputfiles), info = []; return; end% nothing to process
           info = ones(1,length(done));index = find(done==0);N = length(index);
           if N==0, return; end% nothing to process
           shortinputfiles = inputfiles(index);shortoutputfiles = outputfiles(index);
           shortinfo = ones(1,N);
           [parforpath,ID] = openParfor(obj,N);
           parfor i=1:N
               try
                  posepoints = readXMLLandmarks(shortinputfiles{i}.FullName);
                  BatchMapper.saveMatFile(shortoutputfiles{i},posepoints);
               catch
                  shortinfo(i) = 0; 
               end
              BatchMapper.progressParfor; 
           end
           closeParfor(obj,parforpath,ID);
           info(index) = shortinfo;
       end
       function [info,inputfiles] = Obj2MatCleanedFiles(obj,inputpath)
           if nargin<2, inputpath = obj.path31;end
           outputpath = obj.path32;
           [inputfiles,outputfiles,done] = prepIO(obj,inputpath,outputpath,'.obj','.mat');
           if isempty(inputfiles), info = []; return; end% nothing to process
           info = ones(1,length(done));index = find(done==0);N = length(index);
           if N==0, return; end% nothing to process
           shortinputfiles = inputfiles(index);shortoutputfiles = outputfiles(index);
           shortinfo = ones(1,N);
           shape = shape3D;
           [parforpath,ID] = openParfor(obj,N);
           parfor i=1:N
              %i=1;
              try
                  forshape = clone(shape);
                  importWavefront(forshape,shortinputfiles{i}.FileName,shortinputfiles{i}.FullFolder,true);% loading Mevislab obj's so last argument set to true
                  BatchMapper.saveMatFile(shortoutputfiles{i},forshape);
              catch
                  shortinfo(i) = 0;
              end
              BatchMapper.progressParfor;  
           end
           closeParfor(obj,parforpath,ID);
           info(index) = shortinfo;
       end
       function [info,inputfiles] = cleanedMeshObj2shape3D(obj,inputpath)
           [inputfiles,nfiles] = BatchMapper.getFiles(inputpath,'.mat');
           [parforpath,ID] = openParfor(obj,nfiles);
           pathoriginal = obj.path01;
           pathconverted = obj.path11;
           pathcleaned = obj.path32;
           pathposepoints = obj.path21; 
           info = ones(1,nfiles);
           parfor i=1:nfiles
               %i=1;
               try
                   mesh = BatchMapper.loadMatFile(inputfiles{i});
                   shape = convert2Shape3D(mesh);
                   outputfile = clone(inputfiles{i});
                   
%                    % export to original folder 
%                    outputfile.MainFolder = pathoriginal;
%                    outputfile.Extension = '.obj';
%                    BatchMapper.prepFolder(outputfile.FullFolder);
%                    exportWavefront(shape,outputfile.Name,outputfile.FullFolder);
                   
                   % remove texturemap
                   rgb = shape.VertexRGB;
                   shape.TextureMap = [];
                   shape.VertexRGB = rgb;
                   
%                    % export to converted 
%                    outputfile.MainFolder = pathconverted;
%                    outputfile.Extension = '.mat';
%                    BatchMapper.prepFolder(outputfile.FullFolder);
%                    BatchMapper.saveMatFile(outputfile,shape);
                                    
                   % export to cleaned
                   outputfile.MainFolder = pathcleaned;
                   outputfile.Extension = '.mat';
                   BatchMapper.prepFolder(outputfile.FullFolder);
                   BatchMapper.saveMatFile(outputfile,shape);
                   
                   % export pose points
                   if isempty(mesh.PoseLM), BatchMapper.progressParfor; continue; end
                   posepoints = shape3D;
                   posepoints.Vertices = mesh.PoseLM.Vertices';
                   outputfile.MainFolder = pathposepoints;
                   outputfile.Extension = '.mat';
                   BatchMapper.prepFolder(outputfile.FullFolder);
                   if posepoints.nVertices==5
                      BatchMapper.saveMatFile(outputfile,posepoints)
                   else
                      disp('Warning, failed to indicate 5 pose landmarks, results are not saved');
                   end
               catch
                   info(i) = 0;
               end
               BatchMapper.progressParfor;
           end
           closeParfor(obj,parforpath,ID);
       end
   end
   methods (Static = true)% STATIC
       function [out,n] = getRecursiveFolders(path)
                outstr = genpath(path);
                if ~strcmp(computer,'PCWIN64')% windows computer
                    index = strfind(outstr,':');
                else% MAC and/or Ubuntu
                    index = strfind(outstr,';');
                end
                n = length(index);
                out = cell(1,n);
                for i=1:1:n
                     if i==1
                        out{i} = BatchMapper.addBackSlash(outstr(1:index(i)-1));
                     else
                        out{i} = BatchMapper.addBackSlash(outstr(index(i-1)+1:index(i)-1));
                     end
                end
       end
       function [inputfiles,nfiles] = getFiles(inputpath,inputformat)
                [recursiveInput,nRecursiveInput] = BatchMapper.getRecursiveFolders(inputpath);% retrieve recursive folders in input path
                inputfiles = {};counter = 0;l = length(inputpath);% Initialize
                for i=1:1:nRecursiveInput% locate all input files
                    files = dir([recursiveInput{i} '*' inputformat]);% wildcart search on extension
                    for j=1:1:length(files)
                        if strcmp(files(j).name(1:2),'._'), continue; end% these are ghost files to be ommitted.
                        counter = counter+1;
                        inputfiles{counter} = BatchFile; %#ok<*AGROW>% create a Batchfile object, and fill in the properties
                        inputfiles{counter}.FileName = files(j).name;
                        inputfiles{counter}.MainFolder = inputpath;
                        folder = BatchMapper.addBackSlash(files(j).folder);
                        if length(folder)>l
                           inputfiles{counter}.SubFolder = folder(l+1:end);
                        else
                           inputfiles{counter}.SubFolder = [];
                        end
                    end
                end
                nfiles = length(inputfiles);
       end
       function out = addBackSlash(in)
                out = in;
                if ~strcmp(out(end),'/'), out = [out '/'];end
       end
       function prepFolder(in)
           warning off;% due to parallel execution, this can generate warnings
           if ~BatchMapper.myIsFolder(in), mkdir(in);end
           warning on;
       end
       function [folder] = locateFile(path,filename)
                [out,n] = BatchMapper.getRecursiveFolders(path);
                f = 0;
                for i=1:1:n
                    tmp = dir([out{i} filename]);
                    if ~isempty(tmp), f = i; break; end
                end
                if f==0, folder = []; return;end
                folder = out{f};
       end
       function saveMatFile(file,in) %#ok<INUSD>
                save(file.FullName,'in');
       end
       function out = loadMatFile(file)
                in = load(file.FullName);
                loaded = fields(in);
                out = in.(loaded{1});
       end
       function out = mapShape(shape,shapepose,template,templatepose,RM,NRM,display)
            if nargin<7, display = false;end
           % step 0, initialize, making sure the original data is not destroyed
            out = clone(template);out.UserData = [];out.TextureMap = [];
            shape = clone(shape);
            if ~isempty(RM), RM = clone(RM);end
            if ~isempty(NRM), NRM = clone(NRM);end
            cleanManifold(shape);
           % step 1, crude aligment using pose points
            if ~isempty(shapepose)
               [~,~,transform] = procrustes(shapepose.Vertices,templatepose.Vertices,'Scaling',true,'Reflection',false);
               out.Vertices = transform.b*out.Vertices*transform.T + repmat(transform.c(1,:),out.nVertices,1);
            end
           % step2, rigid fine aligment using a rigid mapper
            if ~isempty(RM)
                RM.TargetShape = clone(shape);
                RM.FloatingShape = out;
                RM.Display = display;% make sure there is no visual feedback
                map(RM);
                out = clone(RM.FloatingShape);
            end
           % step 3, non rigid aligment using a non rigid mapper
            if ~isempty(NRM)
                % check on resolution compatibility!
                resFloating = getResolution(out);
                resTarget = getResolution(shape);
                if 1.05*resTarget>=resFloating
                   % The resolution of the target shape is too low, to
                   % allow a complete convergence of smoothing towards 0;
                   % Some final smooting is required, otherwise the outcome
                   % mesh becomes distorted
                   %NRM.TransformNumViscousIterationsEnd = 5;
                   %NRM.TransformNumElasticIterationsEnd = 5;
                   % solution is to upsample the target mesh;
                   NRM.UpSampleTarget = true;
                   NRM.UpSampleMode = 'size';
                   NRM.UpSampleVal = 1.5;
                end
                NRM.TargetShape = clone(shape);
                NRM.FloatingShape = out;
                NRM.Display = display;% maker sure there is no visual feedback
                map(NRM);
                out = clone(NRM.FloatingShape);
                out.UserData.InlierWeights = NRM.InlierWeights;
            end 
       end
       function [shape,shapepose] = reflectShape(shape,shapepose)
           shape.Vertices(:,1) = -1*shape.Vertices(:,1);shape.FlipNormals = true;
           shapepose.Vertices(:,1) = -1*shapepose.Vertices(:,1);
           reindex = [2 1 3 5 4];shapepose.Vertices = shapepose.Vertices(reindex,:);
       end
       function percent = progressParfor(N)
            %PARFOR_PROGRESS Progress monitor (progress bar) that works with parfor.
            %   PARFOR_PROGRESS works by creating a file called parfor_progress.txt in
            %   your working directory, and then keeping track of the parfor loop's
            %   progress within that file. This workaround is necessary because parfor
            %   workers cannot communicate with one another so there is no simple way
            %   to know which iterations have finished and which haven't.
            %
            %   PARFOR_PROGRESS(N) initializes the progress monitor for a set of N
            %   upcoming calculations.
            %
            %   PARFOR_PROGRESS updates the progress inside your parfor loop and
            %   displays an updated progress bar.
            %
            %   PARFOR_PROGRESS(0) deletes parfor_progress.txt and finalizes progress
            %   bar.
            %
            %   To suppress output from any of these functions, just ask for a return
            %   variable from the function calls, like PERCENT = PARFOR_PROGRESS which
            %   returns the percentage of completion.
            %
            %   Example:
            %
            %      N = 100;
            %      parfor_progress(N);
            %      parfor i=1:N
            %         pause(rand); % Replace with real code
            %         parfor_progress;
            %      end
            %      parfor_progress(0);
            %
            %   See also PARFOR.
            % By Jeremy Scheff - jdscheff@gmail.com - http://www.jeremyscheff.com/
            %error(nargchk(0, 1, nargin, 'struct'));

            if nargin < 1,N = -1;end

            percent = 0;
            w = 50; % Width of progress bar

            if N > 0
                f = fopen('parfor_progress.txt', 'w');
                if f<0
                    error('Do you have write permissions for %s?', pwd);
                end
                fprintf(f, '%d\n', N); % Save N at the top of progress.txt
                fclose(f);

                if nargout == 0
                    disp(['  0%[>', repmat(' ', 1, w), ']']);
                end
            elseif N == 0
                delete('parfor_progress.txt');
                percent = 100;

                if nargout == 0
                    disp([repmat(char(8), 1, (w+9)), char(10), '100%[', repmat('=', 1, w+1), ']']);
                end
            else
                if ~exist('parfor_progress.txt', 'file')
                    error('parfor_progress.txt not found. Run PARFOR_PROGRESS(N) before PARFOR_PROGRESS to initialize parfor_progress.txt.');
                end

                f = fopen('parfor_progress.txt', 'a');
                fprintf(f, '1\n');
                fclose(f);

                f = fopen('parfor_progress.txt', 'r');
                progress = fscanf(f, '%d');
                fclose(f);
                percent = (length(progress)-1)/progress(1)*100;

                if nargout == 0
                    perc = sprintf('%3.0f%%', percent); % 4 characters wide, percentage
                    disp([repmat(char(8), 1, (w+9)), char(10), perc, '[', repmat('=', 1, round(percent*w/100)), '>', repmat(' ', 1, w - round(percent*w/100)), ']']);
                end
            end
       end
       function out = alignWithTemplate(shape,template)
           out = clone(shape);
           [~,~,transform] = procrustes(template.Vertices,out.Vertices,'Scaling',true,'Reflection',false);
           out.Vertices = transform.b*out.Vertices*transform.T + repmat(transform.c(1,:),out.nVertices,1);
       end
       function out = myIsFolder(in)
           if verLessThan('matlab', '9.3')
               out = isdir(in);
           else
               out = isfolder(in);
           end
       end       
   end
end