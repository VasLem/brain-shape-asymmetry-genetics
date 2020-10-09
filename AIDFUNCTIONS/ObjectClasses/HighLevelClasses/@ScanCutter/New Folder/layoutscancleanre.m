function varargout = layoutscancleanre(varargin)
% LAYOUTSCANCLEANRE M-file for layoutscancleanre.fig
%      LAYOUTSCANCLEANRE, by itself, creates a new LAYOUTSCANCLEANRE or raises the existing
%      singleton*.
%
%      H = LAYOUTSCANCLEANRE returns the handle to a new LAYOUTSCANCLEANRE or the handle to
%      the existing singleton*.
%
%      LAYOUTSCANCLEANRE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LAYOUTSCANCLEANRE.M with the given input arguments.
%
%      LAYOUTSCANCLEANRE('Property','Value',...) creates a new LAYOUTSCANCLEANRE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before layoutscancleanre_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to layoutscancleanre_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help layoutscancleanre

% Last Modified by GUIDE v2.5 16-Jul-2010 15:22:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @layoutscancleanre_OpeningFcn, ...
                   'gui_OutputFcn',  @layoutscancleanre_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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


% --- Executes just before layoutscancleanre is made visible.
function layoutscancleanre_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to layoutscancleanre (see VARARGIN)

% Choose default command line output for layoutscancleanre
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes layoutscancleanre wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = layoutscancleanre_OutputFcn(hObject, eventdata, handles) 
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


