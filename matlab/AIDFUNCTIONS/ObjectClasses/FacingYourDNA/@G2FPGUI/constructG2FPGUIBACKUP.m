function constructG2FPGUIBACKUP(obj)
% MAIN FIGURE;
obj.Figure = figure(...
                    'Units','characters',...
                    'Position',[3.0000   26.7692   82.4000   32.1538],...
                    'Visible',get(0,'defaultfigureVisible'),...
                    'Color',get(0,'defaultfigureColor'),...
                    'CurrentAxesMode','manual',...
                    'IntegerHandle','off',...
                    'MenuBar','none',...
                    'Name','FacePred V3                     Copyright (c) 2015 P. Claes',...
                    'NumberTitle','off',...
                    'Resize','off',...
                    'PaperPosition',get(0,'defaultfigurePaperPosition'),...
                    'ScreenPixelsPerInchMode','manual',...
                    'ChildrenMode','manual',...
                    'ParentMode','manual',...
                    'HandleVisibility','callback',...
                    'Tag','figure1',...
                    'CloseRequestFcn',@G2FPGUICloseCallback);
setappdata(obj.Figure,'Application',obj);
obj.Handles.LogoAxes = axes(...
                    'Parent',obj.Figure,...
                    'FontUnits',get(0,'defaultaxesFontUnits'),...
                    'Units','characters',...
                    'CameraMode',get(0,'defaultaxesCameraMode'),...
                    'CameraPosition',[0.5 0.5 9.16025403784439],...
                    'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
                    'CameraTarget',[0.5 0.5 0.5],...
                    'CameraTargetMode',get(0,'defaultaxesCameraTargetMode'),...
                    'CameraViewAngle',6.60861036031192,...
                    'CameraViewAngleMode',get(0,'defaultaxesCameraViewAngleMode'),...
                    'Position',[6.8 2.61538461538462 16.6 6.30769230769231],...
                    'ActivePositionProperty','position',...
                    'ActivePositionPropertyMode',get(0,'defaultaxesActivePositionPropertyMode'),...
                    'LooseInset',[14.56 3.55384615384615 10.64 2.42307692307692],...
                    'LooseInsetMode',get(0,'defaultaxesLooseInsetMode'),...
                    'DataSpaceMode',get(0,'defaultaxesDataSpaceMode'),...
                    'PlotBoxAspectRatio',[1 0.987951807228916 0.987951807228916],...
                    'PlotBoxAspectRatioMode',get(0,'defaultaxesPlotBoxAspectRatioMode'),...
                    'ColorSpaceMode',get(0,'defaultaxesColorSpaceMode'),...
                    'ChildContainerMode',get(0,'defaultaxesChildContainerMode'),...
                    'DecorationContainerMode',get(0,'defaultaxesDecorationContainerMode'),...
                    'XRulerMode',get(0,'defaultaxesXRulerMode'),...
                    'XColor',[0.149019607843137 0.149019607843137 0.149019607843137],...
                    'XTick',[0 0.5 1],...
                    'XTickMode',get(0,'defaultaxesXTickMode'),...
                    'XTickLabel',{  '0'; '0.5'; '1' },...
                    'XTickLabelMode',get(0,'defaultaxesXTickLabelMode'),...
                    'XBaselineMode',get(0,'defaultaxesXBaselineMode'),...
                    'YRulerMode',get(0,'defaultaxesYRulerMode'),...
                    'YTick',[0 0.5 1],...
                    'YTickMode',get(0,'defaultaxesYTickMode'),...
                    'YTickLabel',{  '0'; '0.5'; '1' },...
                    'YTickLabelMode',get(0,'defaultaxesYTickLabelMode'),...
                    'YBaselineMode',get(0,'defaultaxesYBaselineMode'),...
                    'ZRulerMode',get(0,'defaultaxesZRulerMode'),...
                    'ZBaselineMode',get(0,'defaultaxesZBaselineMode'),...
                    'AmbientLightSourceMode',get(0,'defaultaxesAmbientLightSourceMode'),...
                    'SortMethod','childorder',...
                    'SortMethodMode',get(0,'defaultaxesSortMethodMode'),...
                    'Tag','axes1',...
                    'ParentMode','manual',...
                    'Visible','off');
obj.Handles.FigStr = uicontrol(...
                    'Parent',obj.Figure,...
                    'FontUnits',get(0,'defaultuicontrolFontUnits'),...
                    'Units','characters',...
                    'String','Peter.Claes@kuleuven.be',...
                    'Style','text',...
                    'Position',[0.8 -0.0769230769230769 28.8 1.76923076923077],...
                    'Children',[],...
                    'ParentMode','manual',...
                    'Tag','text2');
% MENU ITEMS
obj.Handles.FileMenu = uimenu(...
                    'Parent',obj.Figure,...
                    'Label','File',...
                    'ChildrenMode','manual',...
                    'ParentMode','manual',...
                    'Tag','Untitled_1');
obj.Handles.ImportModel = uimenu(...
                    'Parent',obj.Handles.FileMenu,...
                    'Callback',@ImportModelCallback,...
                    'Label','Import Model',...
                    'ParentMode','manual',...
                    'Tag','load');
setappdata(obj.Handles.ImportModel,'Application',obj);
obj.Handles.ImportDNA = uimenu(...
                    'Parent',obj.Handles.FileMenu,...
                    'Callback',@ImportDNACallback,...
                    'Label','Import DNA',...
                    'ParentMode','manual',...
                    'Tag','load');
setappdata(obj.Handles.ImportDNA,'Application',obj);
obj.Handles.Export3D = uimenu(...
                    'Parent',obj.Handles.FileMenu,...
                    'Separator', 'on',...
                    'Callback',@Export3DCallback,...
                    'Label','Export 3D',...
                    'ParentMode','manual',...
                    'Tag','load');
setappdata(obj.Handles.Export3D,'Application',obj);
obj.Handles.Export2D = uimenu(...
                    'Parent',obj.Handles.FileMenu,...
                    'Callback',@Export2DCallback,...
                    'Label','Export 2D',...
                    'ParentMode','manual',...
                    'Tag','load');
setappdata(obj.Handles.Export2D,'Application',obj);
% ACTION PANEL
obj.Handles.ActionPanel = uipanel(...
                    'Parent',obj.Figure,...
                    'FontUnits',get(0,'defaultuipanelFontUnits'),...
                    'Units','characters',...
                    'Title','ACTIONS',...
                    'Position',[3.6 10 36.4 20.6923076923077],...
                    'Visible',get(0,'defaultuipanelVisible'),...
                    'ParentMode','manual',...
                    'Tag','uipanel1');
obj.Handles.PredictButton = uicontrol(...
                    'Parent',obj.Handles.ActionPanel,...
                    'FontUnits',get(0,'defaultuicontrolFontUnits'),...
                    'Units','characters',...
                    'String','FACING DNA',...
                    'Style',get(0,'defaultuicontrolStyle'),...
                    'Position',[6 15.3076923076923 22.6 1.92307692307692],...
                    'Callback',@PredictButtonCallback,...
                    'Children',[],...
                    'ParentMode','manual',...
                    'Tag','pushbutton1');
setappdata(obj.Handles.PredictButton,'Application',obj);
obj.Handles.ClearButton = uicontrol(...
                    'Parent',obj.Handles.ActionPanel,...
                    'FontUnits',get(0,'defaultuicontrolFontUnits'),...
                    'Units','characters',...
                    'String','CLEAR',...
                    'Style',get(0,'defaultuicontrolStyle'),...
                    'Position',[6 12.3076923076923 22.6 1.92307692307692],...
                    'Callback',@ClearButtonCallback,...
                    'Children',[],...
                    'KeyPressFcn',blanks(0),...
                    'ParentMode','manual',...
                    'DeleteFcn',blanks(0),...
                    'ButtonDownFcn',blanks(0),...
                    'Tag','pushbutton2');
setappdata(obj.Handles.ClearButton,'Application',obj);
obj.Handles.SyncButton = uicontrol(...
                    'Parent',obj.Handles.ActionPanel,...
                    'FontUnits',get(0,'defaultuicontrolFontUnits'),...
                    'Units','characters',...
                    'String','SYNC VIEWERS',...
                    'Style',get(0,'defaultuicontrolStyle'),...
                    'Position',[6 9.3076923076923 22.6 1.92307692307692],...
                    'Callback',@SyncButtonCallback,...
                    'Children',[],...
                    'KeyPressFcn',blanks(0),...
                    'ParentMode','manual',...
                    'DeleteFcn',blanks(0),...
                    'ButtonDownFcn',blanks(0),...
                    'Tag','pushbutton2');
setappdata(obj.Handles.SyncButton,'Application',obj);
obj.Handles.AFSlider = uicontrol(...
                    'Parent',obj.Handles.ActionPanel,...
                    'FontUnits',get(0,'defaultuicontrolFontUnits'),...
                    'Units','characters',...
                    'Max',5,...
                    'Min',1,...
                    'String',{  'Slider' },...
                    'Style','slider',...
                    'Value',2,...
                    'Position',[1.8 1.46153846153846 32.4 1.30769230769231],...
                    'BackgroundColor',[0.9 0.9 0.9],...
                    'Callback',@AFSliderCallback,...
                    'Children',[],...
                    'ParentMode','manual',...
                    'Tag','slider1');
setappdata(obj.Handles.AFSlider,'Application',obj);
obj.Handles.StrAF = uicontrol(...
                    'Parent',obj.Handles.ActionPanel,...
                    'FontUnits',get(0,'defaultuicontrolFontUnits'),...
                    'Units','characters',...
                    'String','Amplification Factor',...
                    'Style','text',...
                    'Position',[2 3.38461538461539 31.8 1.30769230769231],...
                    'Children',[],...
                    'Visible',get(0,'defaultuicontrolVisible'),...
                    'ParentMode','manual',...
                    'Tag','text3');
% INFORMATION PANEL
obj.Handles.SettingsPanel = uipanel(...
                    'Parent',obj.Figure,...
                    'FontUnits',get(0,'defaultuipanelFontUnits'),...
                    'Units','characters',...
                    'Title','INFORMATION',...
                    'Position',[42.8 9.92307692307692 36.4 20.6923076923077],...
                    'Visible',get(0,'defaultuipanelVisible'),...
                    'ResizeFcn',blanks(0),...
                    'ParentMode','manual',...
                    'DeleteFcn',blanks(0),...
                    'ButtonDownFcn',blanks(0),...
                    'Tag','uipanel3');
obj.Handles.StrModel = uicontrol(...
                    'Parent',obj.Handles.SettingsPanel,...
                    'FontUnits',get(0,'defaultuicontrolFontUnits'),...
                    'Units','characters',...
                    'String','  ',...
                    'FontSize',8,...
                    'FontWeight','bold',...
                    'ForegroundColor',[0 0 0],...
                    'Style','text',...
                    'Position',[0.8 16 34 2],...
                    'Children',[],...
                    'ParentMode','manual',...
                    'Tag','text4');
obj.Handles.StrSubject = uicontrol(...
                    'Parent',obj.Handles.SettingsPanel,...
                    'FontUnits',get(0,'defaultuicontrolFontUnits'),...
                    'Units','characters',...
                    'String','  ',...
                    'FontSize',8,...
                    'FontWeight','bold',...
                    'ForegroundColor',[0 0 0],...
                    'Style','text',...
                    'Position',[0.8 14 34 2],...
                    'Children',[],...
                    'ParentMode','manual',...
                    'Tag','text5');
obj.Handles.StrProgress = uicontrol(...
                    'Parent',obj.Handles.SettingsPanel,...
                    'FontUnits',get(0,'defaultuicontrolFontUnits'),...
                    'Units','characters',...
                    'String','  ',...
                    'FontSize',8,...
                    'ForegroundColor',[0 0 0],...
                    'Style','text',...
                    'Position',[0.8 1 34 10],...
                    'Children',[],...
                    'ParentMode','manual',...
                    'Tag','text6');                
% VISUALIZATION PANEL
obj.Handles.ViewerPanel = uipanel(...
                    'Parent',obj.Figure,...
                    'FontUnits',get(0,'defaultuipanelFontUnits'),...
                    'Units','characters',...
                    'Title','VIEWERS',...
                    'Position',[29.4 1.923 50.2 7.692],...
                    'Visible',get(0,'defaultuipanelVisible'),...
                    'ResizeFcn',blanks(0),...
                    'ParentMode','manual',...
                    'DeleteFcn',blanks(0),...
                    'ButtonDownFcn',blanks(0),...
                    'Tag','uipanel4');
obj.Handles.BFBox = uicontrol(...
                    'Parent',obj.Handles.ViewerPanel,...
                    'FontUnits',get(0,'defaultuicontrolFontUnits'),...
                    'Units','characters',...
                    'String','Base',...
                    'Style','checkbox',...
                    'Value',1,...
                    'Position',[3 3.84615384615385 13.2 1.76923076923077],...
                    'Callback',@BFBoxCallback,...
                    'Children',[],...
                    'ParentMode','manual',...
                    'Tag','BV');
setappdata(obj.Handles.BFBox,'Application',obj);                
obj.Handles.PFBox = uicontrol(...
                    'Parent',obj.Handles.ViewerPanel,...
                    'FontUnits',get(0,'defaultuicontrolFontUnits'),...
                    'Units','characters',...
                    'String','Prediction',...
                    'Style','checkbox',...
                    'Value',1,...
                    'Position',[3 1.23076923076923 16.8 1.76923076923077],...
                    'Callback',@PFBoxCallback,...
                    'Children',[],...
                    'ParentMode','manual',...
                    'Tag','PV');
setappdata(obj.Handles.PFBox,'Application',obj);
obj.Handles.NDBox = uicontrol(...
                    'Parent',obj.Handles.ViewerPanel,...
                    'FontUnits',get(0,'defaultuicontrolFontUnits'),...
                    'Units','characters',...
                    'String','Feature Prominence',...
                    'Style','checkbox',...
                    'Value',1,...
                    'Position',[21.6 3.76923076923077 26.2 1.76923076923077],...
                    'Callback',@NDBoxCallback,...
                    'Children',[],...
                    'ParentMode','manual',...
                    'Tag','NV');
setappdata(obj.Handles.NDBox,'Application',obj);
obj.Handles.ARBox = uicontrol(...
                    'Parent',obj.Handles.ViewerPanel,...
                    'FontUnits',get(0,'defaultuicontrolFontUnits'),...
                    'Units','characters',...
                    'String','Feature Size',...
                    'Style','checkbox',...
                    'Value',1,...
                    'Position',[21.6 1.23076923076923 26.2 1.76923076923077],...
                    'Callback',@ARBoxCallback,...
                    'Children',[],...
                    'ParentMode','manual',...
                    'Tag','NV');
setappdata(obj.Handles.ARBox,'Application',obj);
obj.Status = 'READY';
end
%% GENERAL CALLBACKS
function G2FPGUICloseCallback(hObject,eventdata) %#ok<*INUSD>
         obj = getappdata(hObject,'Application');
         obj.Status = 'BUSY';
         delete(obj);
         close all;
end
%% MENU CALLBACKS
function ImportModelCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         obj.Handles.StrProgress.String = 'SELECT MODEL FILE';
         [filename, pathname] = uigetfile({'*.mat','PREDICTION MODEL'},'MultiSelect','off');
         if filename == 0, disp('Opertation Cancelled'); return; end
         obj.Status = 'BUSY';
         cd(pathname);
         obj.Handles.StrProgress.String = '...LOADING...';
         drawnow;
         in = load(filename);
         if ~isempty(obj.PM), delete(obj.PM); end
         obj.PM = in.PredModel;
         obj.Handles.StrModel.String = ['Model: ' obj.PM.TrackID];
         obj.Handles.StrProgress.String = 'MODEL LOADED';
         obj.Status = 'READY';
end
function ImportDNACallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         obj.Handles.StrProgress.String = 'SELECT INPUT FILE';
         [filename, pathname] = uigetfile({'*.xlsx','INPUTFILE'},'MultiSelect','off');
         if filename == 0, disp('Opertation Cancelled'); return; end
         obj.Status = 'BUSY';
         cd(pathname);
         disp('  ');
         obj = getappdata(hObject,'Application');
         obj.Subject = clear(obj.Subject);
         parsePredictionInput(obj.Subject,filename);
         obj.Handles.StrSubject.String = ['ID: ' obj.Subject.ID];
         obj.Handles.StrProgress.String = 'INPUT LOADED';
         obj.Status = 'READY';
end
function Export3DCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         obj.Handles.StrProgress.String = 'SELECT EXPORT DIRECTORY';
         obj.Status = 'BUSY';
         obj.Status = 'READY';
end
function Export2DCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         obj.Handles.StrProgress.String = 'SELECT EXPORT DIRECTORY';
         obj.Status = 'BUSY';
         obj.Status = 'READY';
end
%% ACTION CALLBACKS
function PredictButtonCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         obj.Status = 'BUSY';
         str = sprintf('%s \n','1/5 BASEFACE CREATION...  ');
         obj.Handles.StrProgress.String = str;drawnow;
         % BASEFACE
         linkWithBaseCont(obj.Subject,obj.PM.BaseCont);
         obj.Subject.BFType = 'X';
         extractBaseFace(obj.Subject,obj.PM.BaseCont);
         obj.Subject.BaseFace.Axes = obj.BFView.RenderAxes;
         obj.Subject.BaseFace.Visible = true;
         str = sprintf('%s \n','1/5 BASEFACE CREATION...  ','2/5 FACE PREDICTION...    ');
         obj.Handles.StrProgress.String = str;drawnow;
         % PREDICTED FACE
         linkWithSNPCont(obj.Subject,obj.PM.SNPCont);
         extractPredictedFace(obj.Subject,obj.PM.SNPBaseCont,obj.PM.SNPCont);
         str = sprintf('%s \n','1/5 BASEFACE CREATION...  ','2/5 FACE PREDICTION...    ',...
                       '3/5 TEXTURE EXTRACTION... ');
         obj.Handles.StrProgress.String = str;drawnow;
         % TEXTURE EXTRACTION
         obj.Subject.TexMethod = 'BF';
         extractTexture(obj.Subject,obj.PM.TexCont,obj.PM.BaseCont,obj.PM.SNPCont);
         str = sprintf('%s \n','1/5 BASEFACE CREATION...  ','2/5 FACE PREDICTION...    ',...
                       '3/5 TEXTURE EXTRACTION... ','4/5 AMPLIFYING IDENTITY...');
         obj.Handles.StrProgress.String = str;drawnow;
         % AMPLIFYING
         obj.Subject.AF = obj.Handles.AFSlider.Value;
         amplifyIdentity(obj.Subject);
         obj.Subject.AmplFace.Axes = obj.PFView.RenderAxes;
         obj.Subject.AmplFace.Visible = true;
         str = sprintf('%s \n','1/5 BASEFACE CREATION...  ','2/5 FACE PREDICTION...    ',...
                       '3/5 TEXTURE EXTRACTION... ','4/5 AMPLIFYING IDENTITY...',...
                       '5/5 ASSESSING IDENTITY... ');
         obj.Handles.StrProgress.String = str;drawnow;
         % ASSESSING IDENTITY
         assessIdentity(obj.Subject);
         obj.Subject.NDFace.Axes = obj.NDView.RenderAxes;
         set(obj.NDView.RenderAxes,'clim',[-1*max(abs(obj.Subject.NDFace.Value)) max(abs(obj.Subject.NDFace.Value))]);
         colormap(obj.NDView.RenderAxes,'jet');
         obj.Subject.NDFace.Visible = true;
         obj.Subject.ARFace.Axes = obj.ARView.RenderAxes;
         set(obj.ARView.RenderAxes,'clim',[-1*max(abs(obj.Subject.ARFace.Value)) max(abs(obj.Subject.ARFace.Value))]);
         colormap(obj.ARView.RenderAxes,'jet');
         obj.Subject.ARFace.Visible = true;
         str = sprintf('%s \n','1/5 BASEFACE CREATION...  ','2/5 FACE PREDICTION...    ',...
                       '3/5 TEXTURE EXTRACTION... ','4/5 AMPLIFYING IDENTITY...',...
                       '5/5 ASSESSING IDENTITY... ','          DONE            ');
         obj.Handles.StrProgress.String = str;drawnow;
         obj.Status = 'READY';
end
function ClearButtonCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         obj.Status = 'BUSY';
         obj.Handles.StrProgress.String = '...CLEARING...';drawnow;
         obj.Subject = clear(obj.Subject);
         obj.Handles.StrSubject.String = '  ';
         obj.Handles.StrProgress.String = 'INPUT CLEARED';
         obj.Status = 'READY';
end
function SyncButtonCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         obj.Handles.StrProgress.String = '...SYNCING...';drawnow;
         if isempty(obj.BFView), return; end
         if ~isempty(obj.PFView), syncCamera(obj.BFView,obj.PFView);end
         if ~isempty(obj.NDView), syncCamera(obj.BFView,obj.NDView);end
         if ~isempty(obj.ARView), syncCamera(obj.BFView,obj.ARView);end
         obj.Handles.StrProgress.String = 'VIEWERS SYNCED';drawnow;
end
function AFSliderCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         if isempty(obj.Subject.BaseFace), return; end
         if isempty(obj.Subject.PredFace), return; end
         if isempty(obj.Subject.AmplFace), return; end
         dir = obj.Subject.PredFace.Vertices-obj.Subject.BaseFace.Vertices;
         obj.Subject.AmplFace.Vertices = obj.Subject.BaseFace.Vertices + obj.Handles.AFSlider.Value*dir;
end
%% VIEWER CALLBACKS
function BFBoxCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         obj.BFView.Visible = get(hObject,'Value');
end
function PFBoxCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         obj.PFView.Visible = get(hObject,'Value');
end
function NDBoxCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         obj.NDView.Visible = get(hObject,'Value');
end
function ARBoxCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         obj.ARView.Visible = get(hObject,'Value');
end
%%