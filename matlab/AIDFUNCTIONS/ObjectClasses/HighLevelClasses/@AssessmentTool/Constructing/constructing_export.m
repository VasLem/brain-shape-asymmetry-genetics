function varargout = constructing_export(varargin)
% CONSTRUCTING_EXPORT M-file for constructing_export.fig
%      CONSTRUCTING_EXPORT, by itself, creates a new CONSTRUCTING_EXPORT or raises the existing
%      singleton*.
%
%      H = CONSTRUCTING_EXPORT returns the handle to a new CONSTRUCTING_EXPORT or the handle to
%      the existing singleton*.
%
%      CONSTRUCTING_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONSTRUCTING_EXPORT.M with the given input arguments.
%
%      CONSTRUCTING_EXPORT('Property','Value',...) creates a new CONSTRUCTING_EXPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before constructing_export_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to constructing_export_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help constructing_export

% Last Modified by GUIDE v2.5 25-Jun-2010 14:33:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @constructing_export_OpeningFcn, ...
                   'gui_OutputFcn',  @constructing_export_OutputFcn, ...
                   'gui_LayoutFcn',  @constructing_export_LayoutFcn, ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before constructing_export is made visible.
function constructing_export_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to constructing_export (see VARARGIN)

% Choose default command line output for constructing_export
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes constructing_export wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = constructing_export_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function SignificanceEdit_Callback(hObject, eventdata, handles)
% hObject    handle to SignificanceEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SignificanceEdit as text
%        str2double(get(hObject,'String')) returns contents of SignificanceEdit as a double


% --- Executes during object creation, after setting all properties.
function SignificanceEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SignificanceEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in NewButton.
function NewButton_Callback(hObject, eventdata, handles)
% hObject    handle to NewButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in LoadButton.
function LoadButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in SaveButton.
function SaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in RemoveButton.
function RemoveButton_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in RemoveAllButton.
function RemoveAllButton_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveAllButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in SaveAllButton.
function SaveAllButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveAllButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function GlobalThresholdEdit_Callback(hObject, eventdata, handles)
% hObject    handle to GlobalThresholdEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GlobalThresholdEdit as text
%        str2double(get(hObject,'String')) returns contents of GlobalThresholdEdit as a double


% --- Executes during object creation, after setting all properties.
function GlobalThresholdEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GlobalThresholdEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ImportLocalThresholdButton.
function ImportLocalThresholdButton_Callback(hObject, eventdata, handles)
% hObject    handle to ImportLocalThresholdButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in CompensateScaleBox.
function CompensateScaleBox_Callback(hObject, eventdata, handles)
% hObject    handle to CompensateScaleBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CompensateScaleBox


% --- Executes on button press in OutliersOnlyBox.
function OutliersOnlyBox_Callback(hObject, eventdata, handles)
% hObject    handle to OutliersOnlyBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OutliersOnlyBox


% --- Executes on button press in OutliersBinaryBox.
function OutliersBinaryBox_Callback(hObject, eventdata, handles)
% hObject    handle to OutliersBinaryBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OutliersBinaryBox



function RobustnessLevelEdit_Callback(hObject, eventdata, handles)
% hObject    handle to RobustnessLevelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RobustnessLevelEdit as text
%        str2double(get(hObject,'String')) returns contents of RobustnessLevelEdit as a double


% --- Executes during object creation, after setting all properties.
function RobustnessLevelEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RobustnessLevelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function OutliersBinaryEdit_Callback(hObject, eventdata, handles)
% hObject    handle to OutliersBinaryEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OutliersBinaryEdit as text
%        str2double(get(hObject,'String')) returns contents of OutliersBinaryEdit as a double


% --- Executes during object creation, after setting all properties.
function OutliersBinaryEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OutliersBinaryEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ImportScanButton.
function ImportScanButton_Callback(hObject, eventdata, handles)
% hObject    handle to ImportScanButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ImportNormButton.
function ImportNormButton_Callback(hObject, eventdata, handles)
% hObject    handle to ImportNormButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in UpdateButton.
function UpdateButton_Callback(hObject, eventdata, handles)
% hObject    handle to UpdateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in AverageAssessmentButton.
function AverageAssessmentButton_Callback(hObject, eventdata, handles)
% hObject    handle to AverageAssessmentButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in DysmorphogramBox.
function DysmorphogramBox_Callback(hObject, eventdata, handles)
% hObject    handle to DysmorphogramBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DysmorphogramBox


% --- Executes on button press in DistanceMapBox.
function DistanceMapBox_Callback(hObject, eventdata, handles)
% hObject    handle to DistanceMapBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DistanceMapBox


% --- Executes on button press in VectorFieldBox.
function VectorFieldBox_Callback(hObject, eventdata, handles)
% hObject    handle to VectorFieldBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of VectorFieldBox


% --- Executes on button press in ThresholdMapBox.
function ThresholdMapBox_Callback(hObject, eventdata, handles)
% hObject    handle to ThresholdMapBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ThresholdMapBox


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8


% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9




% --- Creates and returns a handle to the GUI figure. 
function h1 = constructing_export_LayoutFcn(policy)
% policy - create a new figure or use a singleton. 'new' or 'reuse'.

persistent hsingleton;
if strcmpi(policy, 'reuse') & ishandle(hsingleton)
    h1 = hsingleton;
    return;
end
load constructing_export.mat


appdata = [];
appdata.GUIDEOptions = struct(...
    'active_h', 216.003540039063, ...
    'taginfo', struct(...
    'figure', 2, ...
    'uipanel', 10, ...
    'listbox', 2, ...
    'uitable', 2, ...
    'pushbutton', 12, ...
    'text', 4, ...
    'edit', 5, ...
    'radiobutton', 3, ...
    'checkbox', 10), ...
    'override', 0, ...
    'release', 13, ...
    'resize', 'none', ...
    'accessibility', 'callback', ...
    'mfile', 1, ...
    'callbacks', 1, ...
    'singleton', 1, ...
    'syscolorfig', 1, ...
    'blocking', 0, ...
    'lastSavedFile', 'C:\MATLAB\Work\ObjectClasses\HighLevelClasses\@AssessmentTool\Constructing\constructing_export.m', ...
    'lastFilename', 'C:\MATLAB\Work\ObjectClasses\HighLevelClasses\@AssessmentTool\Constructing\constructing.fig');
appdata.lastValidTag = 'figure1';
appdata.GUIDELayoutEditor = [];
appdata.initTags = struct(...
    'handle', [], ...
    'tag', 'figure1');

h1 = figure(...
'Units','characters',...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',[0.941176470588235 0.941176470588235 0.941176470588235],...
'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
'IntegerHandle','off',...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'MenuBar','none',...
'Name','constructing',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'PaperSize',[20.98404194812 29.67743169791],...
'PaperType',get(0,'defaultfigurePaperType'),...
'Position',[103.833333333333 14 125.333333333333 47.4375],...
'HandleVisibility','callback',...
'Tag','figure1',...
'UserData',[],...
'Visible','on',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'AssessmentPanel';

h2 = uipanel(...
'Parent',h1,...
'FontSize',10,...
'FontWeight','bold',...
'ForegroundColor',[0 0 0.502],...
'Title','Assessment List',...
'TitlePosition','centertop',...
'Tag','AssessmentPanel',...
'Clipping','on',...
'Position',[0.0127591706539075 0.00810372771474878 0.966507177033493 0.611021069692058],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'AsssesmentTable';

h3 = uitable(...
'Parent',h2,...
'Units','normalized',...
'Data',{  [] []; [] []; [] []; [] [] },...
'Position',[0.0166112956810631 0.106145251396648 0.965116279069767 0.88268156424581],...
'Tag','AsssesmentTable',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'NewButton';

h4 = uicontrol(...
'Parent',h2,...
'Units','normalized',...
'Callback',mat{1},...
'FontWeight','bold',...
'ForegroundColor',[0 0 0.502],...
'Position',[0.0149584487534626 0.0251177394034537 0.119667590027701 0.0613989185417757],...
'String','New',...
'Tag','NewButton',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'LoadButton';

h5 = uicontrol(...
'Parent',h2,...
'Units','normalized',...
'Callback',mat{2},...
'FontWeight','bold',...
'ForegroundColor',[0 0 0.502],...
'Position',[0.146260387811634 0.0251177394034537 0.119667590027701 0.0613989185417757],...
'String','Load',...
'Tag','LoadButton',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'SaveButton';

h6 = uicontrol(...
'Parent',h2,...
'Units','normalized',...
'Callback',mat{3},...
'FontWeight','bold',...
'ForegroundColor',[0 0 0.502],...
'Position',[0.279224376731302 0.0251177394034537 0.119667590027701 0.0613989185417757],...
'String','Save',...
'Tag','SaveButton',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'RemoveButton';

h7 = uicontrol(...
'Parent',h2,...
'Units','normalized',...
'Callback',mat{4},...
'FontWeight','bold',...
'ForegroundColor',[0 0 0.502],...
'Position',[0.731301939058172 0.0251177394034537 0.119667590027701 0.0613989185417757],...
'String','Remove',...
'Tag','RemoveButton',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'RemoveAllButton';

h8 = uicontrol(...
'Parent',h2,...
'Units','normalized',...
'Callback',mat{5},...
'FontWeight','bold',...
'ForegroundColor',[0 0 0.502],...
'Position',[0.860941828254848 0.0279085993371708 0.119667590027701 0.0613989185417757],...
'String','Remove All',...
'Tag','RemoveAllButton',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'SaveAllButton';

h9 = uicontrol(...
'Parent',h2,...
'Units','normalized',...
'Callback',mat{6},...
'FontWeight','bold',...
'ForegroundColor',[0 0 0.502],...
'Position',[0.410526315789474 0.0251177394034537 0.119667590027701 0.0613989185417757],...
'String','Save All',...
'Tag','SaveAllButton',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'Settings';

h10 = uipanel(...
'Parent',h1,...
'FontSize',10,...
'FontWeight','bold',...
'ForegroundColor',[0 0 0.502],...
'Title','Settings',...
'TitlePosition','centertop',...
'Tag','Settings',...
'UserData',[],...
'Clipping','on',...
'Position',[0.0127591706539075 0.623987034035656 0.30622009569378 0.374392220421394],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'uipanel5';

h11 = uibuttongroup(...
'Parent',h10,...
'Title','Threshold Analysis',...
'Tag','uipanel5',...
'Clipping','on',...
'Position',[0.031858407079646 0.0518715001473622 0.929203539823009 0.325375773651635],...
'SelectedObject',[],...
'SelectionChangeFcn',[],...
'OldSelectedObject',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'GlobalRadioButton';

h12 = uicontrol(...
'Parent',h11,...
'Units','normalized',...
'Callback',mat{7},...
'Position',[0.029126213592233 0.53962703962704 0.332038834951456 0.428904428904429],...
'String','Global',...
'Style','radiobutton',...
'Tag','GlobalRadioButton',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'LocalRadioButton';

h13 = uicontrol(...
'Parent',h11,...
'Units','normalized',...
'Callback',mat{8},...
'Position',[0.029126213592233 0.12937062937063 0.332038834951456 0.428904428904429],...
'String','Local',...
'Style','radiobutton',...
'Value',1,...
'Tag','LocalRadioButton',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'GlobalThresholdEdit';

h14 = uicontrol(...
'Parent',h11,...
'Units','normalized',...
'BackgroundColor',[1 1 1],...
'Callback',mat{9},...
'Position',[0.349514563106796 0.632867132867133 0.332038834951456 0.27972027972028],...
'String',blanks(0),...
'Style','edit',...
'CreateFcn',mat{10},...
'Tag','GlobalThresholdEdit');

local_CreateFcn(h14, [], '', appdata);

appdata = [];
appdata.lastValidTag = 'ImportLocalThresholdButton';

h15 = uicontrol(...
'Parent',h11,...
'Units','normalized',...
'Callback',mat{11},...
'FontSize',7,...
'ForegroundColor',[0 0 0.502],...
'Position',[0.337864077669903 0.185314685314686 0.506796116504854 0.372960372960373],...
'String','Import Thresholds',...
'Tag','ImportLocalThresholdButton',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'uipanel6';

h16 = uipanel(...
'Parent',h10,...
'Title','Superimposition',...
'Tag','uipanel6',...
'Clipping','on',...
'Position',[0.031858407079646 0.683760683760684 0.923893805309734 0.311229000884173],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'CompensateScaleBox';

h17 = uicontrol(...
'Parent',h16,...
'Units','normalized',...
'Callback',mat{12},...
'Position',[0.0234146341463415 0.465144230769231 0.913170731707317 0.480769230769231],...
'String','Compensate Scale',...
'Style','checkbox',...
'Tag','CompensateScaleBox',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text3';

h18 = uicontrol(...
'Parent',h16,...
'Units','normalized',...
'Position',[0.298536585365854 0.176682692307693 0.614634146341463 0.288461538461538],...
'String','Robustness Level',...
'Style','text',...
'Tag','text3',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'RobustnessLevelEdit';

h19 = uicontrol(...
'Parent',h16,...
'Units','normalized',...
'BackgroundColor',[1 1 1],...
'Callback',mat{13},...
'Position',[0.0234146341463415 0.157451923076923 0.269268292682927 0.326923076923077],...
'String',blanks(0),...
'Style','edit',...
'CreateFcn',mat{14},...
'Tag','RobustnessLevelEdit');

local_CreateFcn(h19, [], '', appdata);

appdata = [];
appdata.lastValidTag = 'uipanel7';

h20 = uipanel(...
'Parent',h10,...
'Title','Outlier Analysis',...
'Tag','uipanel7',...
'Clipping','on',...
'Position',[0.031858407079646 0.372531682876511 0.907964601769912 0.325375773651635],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'OutliersOnlyBox';

h21 = uicontrol(...
'Parent',h20,...
'Units','normalized',...
'Callback',mat{15},...
'Position',[0.0358208955223881 0.458094144661309 0.93134328358209 0.459242250287026],...
'String','Outliers Only',...
'Style','checkbox',...
'Tag','OutliersOnlyBox',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'OutliersBinaryBox';

h22 = uicontrol(...
'Parent',h20,...
'Units','normalized',...
'Callback',mat{16},...
'Position',[0.0358208955223881 0.109070034443169 0.602985074626866 0.422502870264064],...
'String','Outliers Binary',...
'Style','checkbox',...
'Tag','OutliersBinaryBox',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'OutliersBinaryEdit';

h23 = uicontrol(...
'Parent',h20,...
'Units','normalized',...
'BackgroundColor',[1 1 1],...
'Callback',mat{17},...
'Position',[0.632835820895522 0.182548794489093 0.274626865671642 0.312284730195178],...
'String',blanks(0),...
'Style','edit',...
'CreateFcn',mat{18},...
'Tag','OutliersBinaryEdit');

local_CreateFcn(h23, [], '', appdata);

appdata = [];
appdata.lastValidTag = 'Actions';

h24 = uipanel(...
'Parent',h1,...
'FontSize',10,...
'FontWeight','bold',...
'ForegroundColor',[0 0 0.502],...
'Title','Actions',...
'TitlePosition','centertop',...
'Tag','Actions',...
'Clipping','on',...
'Position',[0.325358851674641 0.625607779578606 0.312599681020734 0.372771474878444],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'uipanel8';

h25 = uipanel(...
'Parent',h24,...
'Title','Single Assessment',...
'Tag','uipanel8',...
'Clipping','on',...
'Position',[0.0467532467532468 0.568047337278107 0.872727272727273 0.421301775147929],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'ImportScanButton';

h26 = uicontrol(...
'Parent',h25,...
'Units','normalized',...
'Callback',mat{19},...
'FontWeight','bold',...
'ForegroundColor',[0 0 0.502],...
'Position',[0.0365482233502538 0.527472527472527 0.505583756345178 0.311073541842773],...
'String','Import Scan',...
'Tag','ImportScanButton',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'ImportNormButton';

h27 = uicontrol(...
'Parent',h25,...
'Units','normalized',...
'Callback',mat{20},...
'FontWeight','bold',...
'ForegroundColor',[0 0 0.502],...
'Position',[0.0426395939086294 0.162299239222316 0.505583756345178 0.311073541842773],...
'String','Import Norm',...
'Tag','ImportNormButton',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'UpdateButton';

h28 = uicontrol(...
'Parent',h25,...
'Units','normalized',...
'Callback',mat{21},...
'FontWeight','bold',...
'ForegroundColor',[0 0 0.502],...
'Position',[0.566497461928934 0.148774302620457 0.395939086294416 0.689771766694844],...
'String','UPDATE',...
'Tag','UpdateButton',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'uipanel9';

h29 = uipanel(...
'Parent',h24,...
'Title','Group Assessments',...
'Tag','uipanel9',...
'Clipping','on',...
'Position',[0.0571428571428571 0.0473372781065089 0.857142857142857 0.501775147928994],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'AverageAssessmentButton';

h30 = uicontrol(...
'Parent',h29,...
'Units','normalized',...
'Callback',mat{22},...
'FontWeight','bold',...
'ForegroundColor',[0 0 0.502],...
'Position',[0.0494845360824742 0.631722260040844 0.896907216494845 0.250510551395507],...
'String','Average Assessment',...
'Tag','AverageAssessmentButton',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'Visualization';

h31 = uipanel(...
'Parent',h1,...
'FontSize',10,...
'FontWeight','bold',...
'ForegroundColor',[0 0 0.502],...
'Title','Visualization',...
'TitlePosition','centertop',...
'Tag','Visualization',...
'Clipping','on',...
'Position',[0.646276595744681 0.623188405797101 0.332446808510638 0.375494071146245],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'DysmorphogramBox';

h32 = uicontrol(...
'Parent',h31,...
'Units','normalized',...
'Callback',mat{23},...
'Position',[0.0682926829268293 0.848758649765738 0.760975609756098 0.117439812096301],...
'String','Dysmorphogram',...
'Style','checkbox',...
'Tag','DysmorphogramBox',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'DistanceMapBox';

h33 = uicontrol(...
'Parent',h31,...
'Units','normalized',...
'Callback',mat{24},...
'Position',[0.0682926829268293 0.731318837669439 0.760975609756098 0.117439812096301],...
'String','Distance Map',...
'Style','checkbox',...
'Tag','DistanceMapBox',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'VectorFieldBox';

h34 = uicontrol(...
'Parent',h31,...
'Units','normalized',...
'Callback',mat{25},...
'Position',[0.0682926829268293 0.613879025573137 0.760975609756098 0.117439812096301],...
'String','Vector Field',...
'Style','checkbox',...
'Tag','VectorFieldBox',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'ThresholdMapBox';

h35 = uicontrol(...
'Parent',h31,...
'Units','normalized',...
'Callback',mat{26},...
'Position',[0.0634146341463415 0.496439213476837 0.760975609756098 0.117439812096301],...
'String','Threshold Map',...
'Style','checkbox',...
'Tag','ThresholdMapBox',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'checkbox8';

h36 = uicontrol(...
'Parent',h31,...
'Units','normalized',...
'Callback',mat{27},...
'Position',[0.0650406504065041 0.381679389312977 0.760162601626017 0.118320610687023],...
'String','Threshold Map',...
'Style','checkbox',...
'Tag','checkbox8',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'checkbox9';

h37 = uicontrol(...
'Parent',h31,...
'Units','normalized',...
'Callback',mat{28},...
'Position',[0.0650406504065041 0.255725190839695 0.760162601626017 0.118320610687023],...
'String','Threshold Map',...
'Style','checkbox',...
'Tag','checkbox9',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );


hsingleton = h1;


% --- Set application data first then calling the CreateFcn. 
function local_CreateFcn(hObject, eventdata, createfcn, appdata)

if ~isempty(appdata)
   names = fieldnames(appdata);
   for i=1:length(names)
       name = char(names(i));
       setappdata(hObject, name, getfield(appdata,name));
   end
end

if ~isempty(createfcn)
   eval(createfcn);
end


% --- Handles default GUIDE GUI creation and callback dispatch
function varargout = gui_mainfcn(gui_State, varargin)

gui_StateFields =  {'gui_Name'
    'gui_Singleton'
    'gui_OpeningFcn'
    'gui_OutputFcn'
    'gui_LayoutFcn'
    'gui_Callback'};
gui_Mfile = '';
for i=1:length(gui_StateFields)
    if ~isfield(gui_State, gui_StateFields{i})
        error('MATLAB:gui_mainfcn:FieldNotFound', 'Could not find field %s in the gui_State struct in GUI M-file %s', gui_StateFields{i}, gui_Mfile);
    elseif isequal(gui_StateFields{i}, 'gui_Name')
        gui_Mfile = [gui_State.(gui_StateFields{i}), '.m'];
    end
end

numargin = length(varargin);

if numargin == 0
    % CONSTRUCTING_EXPORT
    % create the GUI only if we are not in the process of loading it
    % already
    gui_Create = true;
elseif local_isInvokeActiveXCallback(gui_State, varargin{:})
    % CONSTRUCTING_EXPORT(ACTIVEX,...)
    vin{1} = gui_State.gui_Name;
    vin{2} = [get(varargin{1}.Peer, 'Tag'), '_', varargin{end}];
    vin{3} = varargin{1};
    vin{4} = varargin{end-1};
    vin{5} = guidata(varargin{1}.Peer);
    feval(vin{:});
    return;
elseif local_isInvokeHGCallbak(gui_State, varargin{:})
    % CONSTRUCTING_EXPORT('CALLBACK',hObject,eventData,handles,...)
    gui_Create = false;
else
    % CONSTRUCTING_EXPORT(...)
    % create the GUI and hand varargin to the openingfcn
    gui_Create = true;
end

if ~gui_Create
    % In design time, we need to mark all components possibly created in
    % the coming callback evaluation as non-serializable. This way, they
    % will not be brought into GUIDE and not be saved in the figure file
    % when running/saving the GUI from GUIDE.
    designEval = false;
    if (numargin>1 && ishghandle(varargin{2}))
        fig = varargin{2};
        while ~isempty(fig) && ~isa(handle(fig),'figure')
            fig = get(fig,'parent');
        end
        
        designEval = isappdata(0,'CreatingGUIDEFigure') || isprop(fig,'__GUIDEFigure');
    end
        
    if designEval
        beforeChildren = findall(fig);
    end
    
    % evaluate the callback now
    varargin{1} = gui_State.gui_Callback;
    if nargout
        [varargout{1:nargout}] = feval(varargin{:});
    else       
        feval(varargin{:});
    end
    
    % Set serializable of objects created in the above callback to off in
    % design time. Need to check whether figure handle is still valid in
    % case the figure is deleted during the callback dispatching.
    if designEval && ishandle(fig)
        set(setdiff(findall(fig),beforeChildren), 'Serializable','off');
    end
else
    if gui_State.gui_Singleton
        gui_SingletonOpt = 'reuse';
    else
        gui_SingletonOpt = 'new';
    end

    % Check user passing 'visible' P/V pair first so that its value can be
    % used by oepnfig to prevent flickering
    gui_Visible = 'auto';
    gui_VisibleInput = '';
    for index=1:2:length(varargin)
        if length(varargin) == index || ~ischar(varargin{index})
            break;
        end

        % Recognize 'visible' P/V pair
        len1 = min(length('visible'),length(varargin{index}));
        len2 = min(length('off'),length(varargin{index+1}));
        if ischar(varargin{index+1}) && strncmpi(varargin{index},'visible',len1) && len2 > 1
            if strncmpi(varargin{index+1},'off',len2)
                gui_Visible = 'invisible';
                gui_VisibleInput = 'off';
            elseif strncmpi(varargin{index+1},'on',len2)
                gui_Visible = 'visible';
                gui_VisibleInput = 'on';
            end
        end
    end
    
    % Open fig file with stored settings.  Note: This executes all component
    % specific CreateFunctions with an empty HANDLES structure.

    
    % Do feval on layout code in m-file if it exists
    gui_Exported = ~isempty(gui_State.gui_LayoutFcn);
    % this application data is used to indicate the running mode of a GUIDE
    % GUI to distinguish it from the design mode of the GUI in GUIDE. it is
    % only used by actxproxy at this time.   
    setappdata(0,genvarname(['OpenGuiWhenRunning_', gui_State.gui_Name]),1);
    if gui_Exported
        gui_hFigure = feval(gui_State.gui_LayoutFcn, gui_SingletonOpt);

        % make figure invisible here so that the visibility of figure is
        % consistent in OpeningFcn in the exported GUI case
        if isempty(gui_VisibleInput)
            gui_VisibleInput = get(gui_hFigure,'Visible');
        end
        set(gui_hFigure,'Visible','off')

        % openfig (called by local_openfig below) does this for guis without
        % the LayoutFcn. Be sure to do it here so guis show up on screen.
        movegui(gui_hFigure,'onscreen');
    else
        gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt, gui_Visible);
        % If the figure has InGUIInitialization it was not completely created
        % on the last pass.  Delete this handle and try again.
        if isappdata(gui_hFigure, 'InGUIInitialization')
            delete(gui_hFigure);
            gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt, gui_Visible);
        end
    end
    if isappdata(0, genvarname(['OpenGuiWhenRunning_', gui_State.gui_Name]))
        rmappdata(0,genvarname(['OpenGuiWhenRunning_', gui_State.gui_Name]));
    end

    % Set flag to indicate starting GUI initialization
    setappdata(gui_hFigure,'InGUIInitialization',1);

    % Fetch GUIDE Application options
    gui_Options = getappdata(gui_hFigure,'GUIDEOptions');
    % Singleton setting in the GUI M-file takes priority if different
    gui_Options.singleton = gui_State.gui_Singleton;

    if ~isappdata(gui_hFigure,'GUIOnScreen')
        % Adjust background color
        if gui_Options.syscolorfig
            set(gui_hFigure,'Color', get(0,'DefaultUicontrolBackgroundColor'));
        end

        % Generate HANDLES structure and store with GUIDATA. If there is
        % user set GUI data already, keep that also.
        data = guidata(gui_hFigure);
        handles = guihandles(gui_hFigure);
        if ~isempty(handles)
            if isempty(data)
                data = handles;
            else
                names = fieldnames(handles);
                for k=1:length(names)
                    data.(char(names(k)))=handles.(char(names(k)));
                end
            end
        end
        guidata(gui_hFigure, data);
    end

    % Apply input P/V pairs other than 'visible'
    for index=1:2:length(varargin)
        if length(varargin) == index || ~ischar(varargin{index})
            break;
        end

        len1 = min(length('visible'),length(varargin{index}));
        if ~strncmpi(varargin{index},'visible',len1)
            try set(gui_hFigure, varargin{index}, varargin{index+1}), catch break, end
        end
    end

    % If handle visibility is set to 'callback', turn it on until finished
    % with OpeningFcn
    gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
    if strcmp(gui_HandleVisibility, 'callback')
        set(gui_hFigure,'HandleVisibility', 'on');
    end

    feval(gui_State.gui_OpeningFcn, gui_hFigure, [], guidata(gui_hFigure), varargin{:});

    if isscalar(gui_hFigure) && ishandle(gui_hFigure)
        % Handle the default callbacks of predefined toolbar tools in this
        % GUI, if any
        guidemfile('restoreToolbarToolPredefinedCallback',gui_hFigure); 
        
        % Update handle visibility
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);

        % Call openfig again to pick up the saved visibility or apply the
        % one passed in from the P/V pairs
        if ~gui_Exported
            gui_hFigure = local_openfig(gui_State.gui_Name, 'reuse',gui_Visible);
        elseif ~isempty(gui_VisibleInput)
            set(gui_hFigure,'Visible',gui_VisibleInput);
        end
        if strcmpi(get(gui_hFigure, 'Visible'), 'on')
            figure(gui_hFigure);
            
            if gui_Options.singleton
                setappdata(gui_hFigure,'GUIOnScreen', 1);
            end
        end

        % Done with GUI initialization
        if isappdata(gui_hFigure,'InGUIInitialization')
            rmappdata(gui_hFigure,'InGUIInitialization');
        end

        % If handle visibility is set to 'callback', turn it on until
        % finished with OutputFcn
        gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
        if strcmp(gui_HandleVisibility, 'callback')
            set(gui_hFigure,'HandleVisibility', 'on');
        end
        gui_Handles = guidata(gui_hFigure);
    else
        gui_Handles = [];
    end

    if nargout
        [varargout{1:nargout}] = feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    else
        feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    end

    if isscalar(gui_hFigure) && ishandle(gui_hFigure)
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);
    end
end

function gui_hFigure = local_openfig(name, singleton, visible)

% openfig with three arguments was new from R13. Try to call that first, if
% failed, try the old openfig.
if nargin('openfig') == 2
    % OPENFIG did not accept 3rd input argument until R13,
    % toggle default figure visible to prevent the figure
    % from showing up too soon.
    gui_OldDefaultVisible = get(0,'defaultFigureVisible');
    set(0,'defaultFigureVisible','off');
    gui_hFigure = openfig(name, singleton);
    set(0,'defaultFigureVisible',gui_OldDefaultVisible);
else
    gui_hFigure = openfig(name, singleton, visible);
end

function result = local_isInvokeActiveXCallback(gui_State, varargin)

try
    result = ispc && iscom(varargin{1}) ...
             && isequal(varargin{1},gcbo);
catch
    result = false;
end

function result = local_isInvokeHGCallbak(gui_State, varargin)

try
    fhandle = functions(gui_State.gui_Callback);
    result = ~isempty(findstr(gui_State.gui_Name,fhandle.file)) || ...
             (ischar(varargin{1}) ...
             && isequal(ishandle(varargin{2}), 1) ...
             && (~isempty(strfind(varargin{1},[get(varargin{2}, 'Tag'), '_'])) || ...
                ~isempty(strfind(varargin{1}, '_CreateFcn'))) );
catch
    result = false;
end


