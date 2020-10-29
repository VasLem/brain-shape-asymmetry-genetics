function varargout = constructing(varargin)
% CONSTRUCTING M-file for constructing.fig
%      CONSTRUCTING, by itself, creates a new CONSTRUCTING or raises the existing
%      singleton*.
%
%      H = CONSTRUCTING returns the handle to a new CONSTRUCTING or the handle to
%      the existing singleton*.
%
%      CONSTRUCTING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONSTRUCTING.M with the given input arguments.
%
%      CONSTRUCTING('Property','Value',...) creates a new CONSTRUCTING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before constructing_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to constructing_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help constructing

% Last Modified by GUIDE v2.5 25-Jun-2010 14:31:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @constructing_OpeningFcn, ...
                   'gui_OutputFcn',  @constructing_OutputFcn, ...
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


% --- Executes just before constructing is made visible.
function constructing_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to constructing (see VARARGIN)

% Choose default command line output for constructing
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes constructing wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = constructing_OutputFcn(hObject, eventdata, handles) 
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


