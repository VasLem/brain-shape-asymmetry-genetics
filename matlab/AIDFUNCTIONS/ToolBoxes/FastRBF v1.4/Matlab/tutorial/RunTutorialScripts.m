% FILENAME: RunTutorialScripts.m
function x = runtutorialscripts(varargin)
% Runs the tutorial scripts
% Usage: runtutorialscripts([startsection, [endsection]])
%        where start and end are chapter.section numbers to be processed
%        e.g., runtutorialscripts(6.3) will run the scripts for section 6.3
%        or
%        runtutorialscripts(5.3, 7.1) will run all the scripts for sections
%        5.3 to 7.1 inclusive
%        or
%        runtutorialscripts(sections) where sections is an array of target
%        sections
%        Current chapter.sections are 5.2:4, 6.1:4, 7.1:3, 8.1:2

% process function arguments
maxchapter = 10;
maxsection = 9;
if nargin == 0
    sections = [1:0.1:10.9];
elseif nargin == 1
    sections = varargin{1};
    [r, n] = size(sections);
    if n == 1
        % single section has been specified
        sections = [sections(1,1)];
    else
        % set of sections has been specified
    end
else
    sections = [varargin{1}:0.1:varargin{2}];
    %disp(eval('sections'))
end

% give access to figure_by_name function
addpath(eval('pwd'));

% Test tutorial examples

%% CHAPTER 5
% section 5.2
if find(sections==5.2) > 0
    cd ImportData;
    ImportAndExportFormatTextData;
    cd ..;
end

% section 5.3
if find(sections==5.3) > 0
    cd ImportData;
    ImportingAndExporting3Ddata;
    cd ..;
end    

% section 5.4
if find(sections==5.4) > 0
    cd ImportData;
    ImportAndExportSurfaceData;
    cd ..;
end

%% CHAPTER 6
% section 6.1
if find(sections==6.1) > 0
    cd FitAndEval;
    FittingAndEvaluating2Ddata;
    cd ..;
end

% section 6.2
if find(sections==6.2) > 0
    cd FitAndEval;
    FittingAndEvaluating3Ddata;
    cd ..;
end

% section 6.3
if find(sections==6.3) > 0
    cd FitAndEval;
    Anti_Aliasing;
    cd ..;
end

% section 6.4
if find(sections==6.4) > 0
    cd FitAndEval;
    NoiseReduction;
    cd ..;
end

%% CHAPTER 7
% section 7.1
if find(sections==7.1) > 0
    cd SurfaceFit;
    FittingSurfacesToMeshData;
    cd ..;
end

% section 7.2
if find(sections==7.2) > 0
    cd ColourSurf;
    ColourSurfaceFitting;
    cd ..;
end

% section 7.3
if find(sections==7.3) > 0
    cd Smoothing;
    SurfacingNoisyData;
    cd ..;
end

%% CHAPTER 8
% section 8.1
if find(sections==8.1) > 0
    cd Isosurf;
    Isosurfacing;
    cd ..;
end

% section 8.2
if find(sections==8.2) > 0
    cd Simplify;
    MeshSimplifying;
    cd ..;
end

rmpath(eval('pwd'));
