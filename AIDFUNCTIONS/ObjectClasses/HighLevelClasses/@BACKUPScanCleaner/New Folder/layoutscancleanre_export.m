function varargout = layoutscancleanre_export(varargin)
% LAYOUTSCANCLEANRE_EXPORT M-file for layoutscancleanre_export.fig
%      LAYOUTSCANCLEANRE_EXPORT, by itself, creates a new LAYOUTSCANCLEANRE_EXPORT or raises the existing
%      singleton*.
%
%      H = LAYOUTSCANCLEANRE_EXPORT returns the handle to a new LAYOUTSCANCLEANRE_EXPORT or the handle to
%      the existing singleton*.
%
%      LAYOUTSCANCLEANRE_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LAYOUTSCANCLEANRE_EXPORT.M with the given input arguments.
%
%      LAYOUTSCANCLEANRE_EXPORT('Property','Value',...) creates a new LAYOUTSCANCLEANRE_EXPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before layoutscancleanre_export_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to layoutscancleanre_export_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help layoutscancleanre_export

% Last Modified by GUIDE v2.5 16-Jul-2010 15:24:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @layoutscancleanre_export_OpeningFcn, ...
                   'gui_OutputFcn',  @layoutscancleanre_export_OutputFcn, ...
                   'gui_LayoutFcn',  @layoutscancleanre_export_LayoutFcn, ...
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


% --- Executes just before layoutscancleanre_export is made visible.
function layoutscancleanre_export_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to layoutscancleanre_export (see VARARGIN)

% Choose default command line output for layoutscancleanre_export
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes layoutscancleanre_export wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = layoutscancleanre_export_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function InputPathEdit_Callback(hObject, eventdata, handles)
% hObject    handle to InputPathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of InputPathEdit as text
%        str2double(get(hObject,'String')) returns contents of InputPathEdit as a double


% --- Executes during object creation, after setting all properties.
function InputPathEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InputPathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in InputPathButton.
function InputPathButton_Callback(hObject, eventdata, handles)
% hObject    handle to InputPathButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in IncludeSubFoldersBox.
function IncludeSubFoldersBox_Callback(hObject, eventdata, handles)
% hObject    handle to IncludeSubFoldersBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of IncludeSubFoldersBox


% --- Executes on selection change in InputFormatMenu.
function InputFormatMenu_Callback(hObject, eventdata, handles)
% hObject    handle to InputFormatMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns InputFormatMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from InputFormatMenu


% --- Executes during object creation, after setting all properties.
function InputFormatMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InputFormatMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function OutputPathEdit_Callback(hObject, eventdata, handles)
% hObject    handle to OutputPathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OutputPathEdit as text
%        str2double(get(hObject,'String')) returns contents of OutputPathEdit as a double


% --- Executes during object creation, after setting all properties.
function OutputPathEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OutputPathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in OutputPathButton.
function OutputPathButton_Callback(hObject, eventdata, handles)
% hObject    handle to OutputPathButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in OutputFormatMenu.
function OutputFormatMenu_Callback(hObject, eventdata, handles)
% hObject    handle to OutputFormatMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns OutputFormatMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from OutputFormatMenu


% --- Executes during object creation, after setting all properties.
function OutputFormatMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OutputFormatMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CenterVerticesBox.
function CenterVerticesBox_Callback(hObject, eventdata, handles)
% hObject    handle to CenterVerticesBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CenterVerticesBox


% --- Executes on button press in DeleteBorderBox.
function DeleteBorderBox_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteBorderBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DeleteBorderBox



function DeleteBorderEdit_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteBorderEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DeleteBorderEdit as text
%        str2double(get(hObject,'String')) returns contents of DeleteBorderEdit as a double


% --- Executes during object creation, after setting all properties.
function DeleteBorderEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DeleteBorderEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PreViewBox.
function PreViewBox_Callback(hObject, eventdata, handles)
% hObject    handle to PreViewBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PreViewBox


% --- Executes on button press in DownSampleBox.
function DownSampleBox_Callback(hObject, eventdata, handles)
% hObject    handle to DownSampleBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DownSampleBox



function DownSampleEdit_Callback(hObject, eventdata, handles)
% hObject    handle to DownSampleEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DownSampleEdit as text
%        str2double(get(hObject,'String')) returns contents of DownSampleEdit as a double


% --- Executes during object creation, after setting all properties.
function DownSampleEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DownSampleEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in IndicatePoseLMBox.
function IndicatePoseLMBox_Callback(hObject, eventdata, handles)
% hObject    handle to IndicatePoseLMBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of IndicatePoseLMBox


% --- Executes on button press in CropBox.
function CropBox_Callback(hObject, eventdata, handles)
% hObject    handle to CropBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CropBox



function CropEdit_Callback(hObject, eventdata, handles)
% hObject    handle to CropEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CropEdit as text
%        str2double(get(hObject,'String')) returns contents of CropEdit as a double


% --- Executes during object creation, after setting all properties.
function CropEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CropEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in DownSampleTextureMapBox.
function DownSampleTextureMapBox_Callback(hObject, eventdata, handles)
% hObject    handle to DownSampleTextureMapBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DownSampleTextureMapBox


% --- Executes on button press in PostViewBox.
function PostViewBox_Callback(hObject, eventdata, handles)
% hObject    handle to PostViewBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PostViewBox


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in OverwriteBox.
function OverwriteBox_Callback(hObject, eventdata, handles)
% hObject    handle to OverwriteBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OverwriteBox




% --- Creates and returns a handle to the GUI figure. 
function h1 = layoutscancleanre_export_LayoutFcn(policy)
% policy - create a new figure or use a singleton. 'new' or 'reuse'.

persistent hsingleton;
if strcmpi(policy, 'reuse') & ishandle(hsingleton)
    h1 = hsingleton;
    return;
end
load layoutscancleanre_export.mat


appdata = [];
appdata.GUIDEOptions = struct(...
    'active_h', [], ...
    'taginfo', struct(...
    'figure', 2, ...
    'uipanel', 3, ...
    'text', 5, ...
    'edit', 6, ...
    'pushbutton', 4, ...
    'checkbox', 11, ...
    'popupmenu', 3), ...
    'override', 0, ...
    'release', 13, ...
    'resize', 'none', ...
    'accessibility', 'callback', ...
    'mfile', 1, ...
    'callbacks', 1, ...
    'singleton', 1, ...
    'syscolorfig', 1, ...
    'blocking', 0, ...
    'lastSavedFile', 'C:\MATLAB\Work\ObjectClasses\HighLevelClasses\@ScanCleaner\New Folder\layoutscancleanre_export.m', ...
    'lastFilename', 'C:\MATLAB\Work\ObjectClasses\HighLevelClasses\@ScanCleaner\New Folder\layoutscancleanre.fig');
appdata.lastValidTag = 'figure1';
appdata.GUIDELayoutEditor = [];
appdata.initTags = struct(...
    'handle', [], ...
    'tag', 'figure1');

h1 = figure(...
'Units','normalized',...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',[0.941176470588235 0.941176470588235 0.941176470588235],...
'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
'IntegerHandle','off',...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'MenuBar','none',...
'Name','layoutscancleanre',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'PaperSize',[20.98404194812 29.67743169791],...
'PaperType',get(0,'defaultfigurePaperType'),...
'Position',[0.2703125 0.101666666666667 0.2171875 0.564166666666667],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1',...
'UserData',[],...
'Visible','on',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'InputOutputPanel';

h2 = uipanel(...
'Parent',h1,...
'FontAngle','oblique',...
'FontSize',10,...
'FontWeight','bold',...
'ForegroundColor',[0.392156862745098 0.474509803921569 0.635294117647059],...
'ShadowColor',[0.305882352941176 0.396078431372549 0.580392156862745],...
'Title','Input & Output',...
'TitlePosition','centertop',...
'Tag','InputOutputPanel',...
'Clipping','on',...
'Position',[0.0239808153477218 0.58493353028065 0.935251798561151 0.389955686853767],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text1';

h3 = uicontrol(...
'Parent',h2,...
'Units','normalized',...
'FontSize',10,...
'FontWeight','bold',...
'Position',[0.0181347150259067 0.865284974093264 0.261658031088083 0.0984455958549225],...
'String','Input Directory:',...
'Style','text',...
'Tag','text1',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'InputPathEdit';

h4 = uicontrol(...
'Parent',h2,...
'Units','normalized',...
'BackgroundColor',[1 1 1],...
'Callback',mat{1},...
'Position',[0.0233160621761658 0.765422438405416 0.738341968911917 0.108808290155441],...
'String',blanks(0),...
'Style','edit',...
'CreateFcn',mat{2},...
'Tag','InputPathEdit');

local_CreateFcn(h4, [], '', appdata);

appdata = [];
appdata.lastValidTag = 'InputPathButton';

h5 = uicontrol(...
'Parent',h2,...
'Units','normalized',...
'Callback',mat{3},...
'FontAngle','oblique',...
'FontWeight','bold',...
'ForegroundColor',[0.392156862745098 0.474509803921569 0.635294117647059],...
'Position',[0.772020725388601 0.763265306122454 0.207253886010363 0.126530612244898],...
'String','Browse',...
'Tag','InputPathButton',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'IncludeSubFoldersBox';

h6 = uicontrol(...
'Parent',h2,...
'Units','normalized',...
'Callback',mat{4},...
'FontSize',10,...
'Position',[0.0310880829015544 0.636734693877552 0.458549222797927 0.13469387755102],...
'String','Include SubFolders',...
'Style','checkbox',...
'Tag','IncludeSubFoldersBox',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text2';

h7 = uicontrol(...
'Parent',h2,...
'Units','normalized',...
'FontSize',10,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[0.0259067357512953 0.534693877551021 0.261658031088083 0.0979591836734694],...
'String','Input Format',...
'Style','text',...
'Tag','text2',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'InputFormatMenu';

h8 = uicontrol(...
'Parent',h2,...
'Units','normalized',...
'BackgroundColor',[1 1 1],...
'Callback',mat{5},...
'Position',[0.297927461139896 0.514285714285714 0.580310880829016 0.118367346938776],...
'String',blanks(0),...
'Style','popupmenu',...
'Value',1,...
'CreateFcn',mat{6},...
'Tag','InputFormatMenu');

local_CreateFcn(h8, [], '', appdata);

appdata = [];
appdata.lastValidTag = 'text3';

h9 = uicontrol(...
'Parent',h2,...
'Units','normalized',...
'FontSize',10,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[0.0233160621761658 0.371428571428573 0.411917098445596 0.122448979591837],...
'String','Output Directory:',...
'Style','text',...
'Tag','text3',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'OutputPathEdit';

h10 = uicontrol(...
'Parent',h2,...
'Units','normalized',...
'BackgroundColor',[1 1 1],...
'Callback',mat{7},...
'Position',[0.0233160621761658 0.281632653061225 0.738341968911917 0.110204081632653],...
'String',blanks(0),...
'Style','edit',...
'CreateFcn',mat{8},...
'Tag','OutputPathEdit');

local_CreateFcn(h10, [], '', appdata);

appdata = [];
appdata.lastValidTag = 'OutputPathButton';

h11 = uicontrol(...
'Parent',h2,...
'Units','normalized',...
'Callback',mat{9},...
'FontAngle','oblique',...
'FontWeight','bold',...
'ForegroundColor',[0.392156862745098 0.474509803921569 0.635294117647059],...
'Position',[0.769430051813472 0.269387755102042 0.207253886010363 0.126530612244898],...
'String','Browse',...
'Tag','OutputPathButton',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text4';

h12 = uicontrol(...
'Parent',h2,...
'Units','normalized',...
'FontSize',10,...
'FontWeight','bold',...
'Position',[0.0181347150259067 0.0408163265306125 0.261658031088083 0.0979591836734694],...
'String','Output Format',...
'Style','text',...
'Tag','text4',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'OutputFormatMenu';

h13 = uicontrol(...
'Parent',h2,...
'Units','normalized',...
'BackgroundColor',[1 1 1],...
'Callback',mat{10},...
'Position',[0.295336787564767 0.0163265306122451 0.580310880829016 0.118367346938776],...
'String',blanks(0),...
'Style','popupmenu',...
'Value',1,...
'CreateFcn',mat{11},...
'Tag','OutputFormatMenu');

local_CreateFcn(h13, [], '', appdata);

appdata = [];
appdata.lastValidTag = 'OverwriteBox';

h14 = uicontrol(...
'Parent',h2,...
'Units','normalized',...
'Callback',mat{12},...
'FontSize',10,...
'Position',[0.0284974093264249 0.138775510204082 0.458549222797927 0.13469387755102],...
'String','Overwrite',...
'Style','checkbox',...
'Tag','OverwriteBox',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'uipanel2';

h15 = uipanel(...
'Parent',h1,...
'FontAngle','oblique',...
'FontSize',10,...
'FontWeight','bold',...
'ForegroundColor',[0.392156862745098 0.474509803921569 0.635294117647059],...
'Title','Actions',...
'TitlePosition','centertop',...
'Tag','uipanel2',...
'Clipping','on',...
'Position',[0.0311750599520384 0.0221565731166913 0.928057553956834 0.556868537666174],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'CenterVerticesBox';

h16 = uicontrol(...
'Parent',h15,...
'Units','normalized',...
'Callback',mat{13},...
'FontSize',10,...
'Position',[0.0365535248041776 0.902234636871509 0.399477806788512 0.0754189944134078],...
'String','Center Vertices',...
'Style','checkbox',...
'Tag','CenterVerticesBox',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'DeleteBorderBox';

h17 = uicontrol(...
'Parent',h15,...
'Units','normalized',...
'Callback',mat{14},...
'FontSize',10,...
'Position',[0.0365535248041776 0.731843575418994 0.399477806788512 0.0754189944134078],...
'String','Detete Border',...
'Style','checkbox',...
'Tag','DeleteBorderBox',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'DeleteBorderEdit';

h18 = uicontrol(...
'Parent',h15,...
'Units','normalized',...
'BackgroundColor',[1 1 1],...
'Callback',mat{15},...
'Position',[0.352480417754569 0.740223463687151 0.195822454308094 0.0726256983240224],...
'String',blanks(0),...
'Style','edit',...
'CreateFcn',mat{16},...
'Tag','DeleteBorderEdit');

local_CreateFcn(h18, [], '', appdata);

appdata = [];
appdata.lastValidTag = 'PreViewBox';

h19 = uicontrol(...
'Parent',h15,...
'Units','normalized',...
'Callback',mat{17},...
'FontSize',10,...
'Position',[0.0365535248041776 0.82122905027933 0.399477806788512 0.0754189944134078],...
'String','Pre View',...
'Style','checkbox',...
'Tag','PreViewBox',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'DownSampleBox';

h20 = uicontrol(...
'Parent',h15,...
'Units','normalized',...
'Callback',mat{18},...
'FontSize',10,...
'Position',[0.0391644908616188 0.639664804469274 0.399477806788512 0.0754189944134078],...
'String','Down Sample',...
'Style','checkbox',...
'Tag','DownSampleBox',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'DownSampleEdit';

h21 = uicontrol(...
'Parent',h15,...
'Units','normalized',...
'BackgroundColor',[1 1 1],...
'Callback',mat{19},...
'Position',[0.352480417754569 0.648044692737431 0.195822454308094 0.0726256983240224],...
'String',blanks(0),...
'Style','edit',...
'CreateFcn',mat{20},...
'Tag','DownSampleEdit');

local_CreateFcn(h21, [], '', appdata);

appdata = [];
appdata.lastValidTag = 'IndicatePoseLMBox';

h22 = uicontrol(...
'Parent',h15,...
'Units','normalized',...
'Callback',mat{21},...
'FontSize',10,...
'Position',[0.0391644908616188 0.541899441340782 0.540469973890339 0.0865921787709497],...
'String','Indicate Pose Landmarks',...
'Style','checkbox',...
'Tag','IndicatePoseLMBox',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'CropBox';

h23 = uicontrol(...
'Parent',h15,...
'Units','normalized',...
'Callback',mat{22},...
'FontSize',10,...
'Position',[0.0391644908616188 0.455307262569832 0.399477806788512 0.0754189944134078],...
'String','Crop',...
'Style','checkbox',...
'Tag','CropBox',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'CropEdit';

h24 = uicontrol(...
'Parent',h15,...
'Units','normalized',...
'BackgroundColor',[1 1 1],...
'Callback',mat{23},...
'Position',[0.35509138381201 0.463687150837989 0.195822454308094 0.0726256983240224],...
'String',blanks(0),...
'Style','edit',...
'CreateFcn',mat{24},...
'Tag','CropEdit');

local_CreateFcn(h24, [], '', appdata);

appdata = [];
appdata.lastValidTag = 'DownSampleTextureMapBox';

h25 = uicontrol(...
'Parent',h15,...
'Units','normalized',...
'Callback',mat{25},...
'FontSize',10,...
'Position',[0.0391644908616188 0.360335195530726 0.577023498694517 0.0865921787709497],...
'String','Down Sample TextureMap',...
'Style','checkbox',...
'Tag','DownSampleTextureMapBox',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'PostViewBox';

h26 = uicontrol(...
'Parent',h15,...
'Units','normalized',...
'Callback',mat{26},...
'FontSize',10,...
'Position',[0.0417754569190601 0.273743016759777 0.399477806788512 0.0754189944134078],...
'String','Post View',...
'Style','checkbox',...
'Tag','PostViewBox',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'pushbutton3';

h27 = uicontrol(...
'Parent',h15,...
'Units','normalized',...
'Callback',mat{27},...
'FontAngle','oblique',...
'FontSize',16,...
'FontWeight','bold',...
'ForegroundColor',[0.392156862745098 0.474509803921569 0.635294117647059],...
'Position',[0.0391644908616188 0.0251396648044693 0.919060052219321 0.201117318435754],...
'String','Start Batch Process',...
'Tag','pushbutton3',...
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
    % LAYOUTSCANCLEANRE_EXPORT
    % create the GUI only if we are not in the process of loading it
    % already
    gui_Create = true;
elseif local_isInvokeActiveXCallback(gui_State, varargin{:})
    % LAYOUTSCANCLEANRE_EXPORT(ACTIVEX,...)
    vin{1} = gui_State.gui_Name;
    vin{2} = [get(varargin{1}.Peer, 'Tag'), '_', varargin{end}];
    vin{3} = varargin{1};
    vin{4} = varargin{end-1};
    vin{5} = guidata(varargin{1}.Peer);
    feval(vin{:});
    return;
elseif local_isInvokeHGCallbak(gui_State, varargin{:})
    % LAYOUTSCANCLEANRE_EXPORT('CALLBACK',hObject,eventData,handles,...)
    gui_Create = false;
else
    % LAYOUTSCANCLEANRE_EXPORT(...)
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


