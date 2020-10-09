function constructScanMapper(obj)
% GENERATING THE MAIN FIGURE PANEL
  obj.Handles.Panel = figure(...
'Units','normalized',...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',[0.941176470588235 0.941176470588235 0.941176470588235],...
'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
'IntegerHandle','off',...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'MenuBar','none',...
'Name',obj.Tag,...
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
'CloseRequestFcn',@PanelCloseCallback);
setappdata(obj.Handles.Panel,'Application',obj);
obj.Handles.InputOutputPanel = uipanel(...
'Parent',obj.Handles.Panel,...
'FontAngle','oblique',...
'FontSize',10,...
'FontWeight','bold',...
'ForegroundColor',[0 0 0],...
'ShadowColor',[0.305882352941176 0.396078431372549 0.580392156862745],...
'Title','Input & Output',...
'TitlePosition','centertop',...
'Tag','InputOutputPanel',...
'Clipping','on',...
'Position',[0.0239808153477218 0.58493353028065 0.935251798561151 0.389955686853767],...
'CreateFcn', {} );
 uicontrol(...
'Parent',obj.Handles.InputOutputPanel,...
'Units','normalized',...
'FontSize',10,...
'FontWeight','bold',...
'Position',[0.0181347150259067 0.865284974093264 0.261658031088083 0.0984455958549225],...
'String','Input Directory:',...
'Style','text',...
'Tag','text1',...
'CreateFcn', {} );
obj.Handles.InputPathEdit = uicontrol(...
'Parent',obj.Handles.InputOutputPanel,...
'Units','normalized',...
'BackgroundColor',[1 1 1],...
'Callback',{},...
'Position',[0.0233160621761658 0.765422438405416 0.738341968911917 0.108808290155441],...
'String',obj.B.InputFolder,...
'Style','edit',...
'CreateFcn',{},...
'Tag','InputPathEdit');
obj.Handles.InputPathButton = uicontrol(...
'Parent',obj.Handles.InputOutputPanel,...
'Units','normalized',...
'Callback',@InputPathButtonCallback,...
'FontAngle','oblique',...
'FontWeight','bold',...
'ForegroundColor',[0 0 0],...
'Position',[0.772020725388601 0.763265306122454 0.207253886010363 0.126530612244898],...
'String','Browse',...
'Tag','InputPathButton',...
'CreateFcn', {} );
 setappdata(obj.Handles.InputPathButton,'Application',obj);
 obj.Handles.IncludeSubFoldersBox = uicontrol(...
'Parent',obj.Handles.InputOutputPanel,...
'Units','normalized',...
'Callback',{},...
'FontSize',10,...
'Position',[0.0310880829015544 0.63 0.458549222797927 0.13469387755102],...
'String','Include SubFolders',...
'Style','checkbox',...
'Tag','IncludeSubFoldersBox',...
'CreateFcn', {} );
uicontrol(...
'Parent',obj.Handles.InputOutputPanel,...
'Units','normalized',...
'FontSize',10,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[0.0259067357512953 0.534693877551021 0.261658031088083 0.0979591836734694],...
'String','Input Format:',...
'Style','text',...
'Tag','text2',...
'CreateFcn', {} );
obj.Handles.InputFormatMenu = uicontrol(...
'Parent',obj.Handles.InputOutputPanel,...
'Units','normalized',...
'BackgroundColor',[1 1 1],...
'Callback',{},...
'Position',[0.297927461139896 0.514285714285714 0.580310880829016 0.118367346938776],...
'String',obj.B.ImportFormats,...
'Style','popupmenu',...
'Value',1,...
'CreateFcn',{},...
'Tag','InputFormatMenu');
uicontrol(...
'Parent',obj.Handles.InputOutputPanel,...
'Units','normalized',...
'FontSize',10,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[0.0233160621761658 0.371428571428573 0.411917098445596 0.122448979591837],...
'String','Output Directory:',...
'Style','text',...
'Tag','text3',...
'CreateFcn', {} );
obj.Handles.OutputPathEdit = uicontrol(...
'Parent',obj.Handles.InputOutputPanel,...
'Units','normalized',...
'BackgroundColor',[1 1 1],...
'Callback',{},...
'Position',[0.0233160621761658 0.281632653061225 0.738341968911917 0.110204081632653],...
'String',blanks(0),...
'Style','edit',...
'CreateFcn',{},...
'Tag','OutputPathEdit');
obj.Handles.OutputPathButton = uicontrol(...
'Parent',obj.Handles.InputOutputPanel,...
'Units','normalized',...
'Callback',@OutputPathButtonCallback,...
'FontAngle','oblique',...
'FontWeight','bold',...
'ForegroundColor',[0 0 0],...
'Position',[0.769430051813472 0.269387755102042 0.207253886010363 0.126530612244898],...
'String','Browse',...
'Tag','OutputPathButton',...
'CreateFcn', {} );
 setappdata(obj.Handles.OutputPathButton,'Application',obj);
 uicontrol(...
'Parent',obj.Handles.InputOutputPanel,...
'Units','normalized',...
'FontSize',10,...
'FontWeight','bold',...
'Position',[0.0181347150259067 0.0408163265306125 0.261658031088083 0.0979591836734694],...
'String','Output Format',...
'Style','text',...
'Tag','text4',...
'CreateFcn', {} );
obj.Handles.OutputFormatMenu = uicontrol(...
'Parent',obj.Handles.InputOutputPanel,...
'Units','normalized',...
'BackgroundColor',[1 1 1],...
'Callback',{},...
'Position',[0.295336787564767 0.0163265306122451 0.580310880829016 0.118367346938776],...
'String',obj.B.ExportFormats,...
'Style','popupmenu',...
'Value',1,...
'CreateFcn',{},...
'Tag','OutputFormatMenu');
obj.Handles.OverwriteBox = uicontrol(...
'Parent',obj.Handles.InputOutputPanel,...
'Units','normalized',...
'Callback',{},...
'FontSize',10,...
'Position',[0.0284974093264249 0.138775510204082 0.458549222797927 0.13469387755102],...
'String','Overwrite',...
'Style','checkbox',...
'Tag','OverwriteBox',...
'CreateFcn', {} );
% ACTIONS
obj.Handles.ActionPanel = uipanel(...
'Parent',obj.Handles.Panel,...
'FontAngle','oblique',...
'FontSize',10,...
'FontWeight','bold',...
'ForegroundColor',[0 0 0],...
'Title','Actions',...
'TitlePosition','centertop',...
'Tag','uipanel2',...
'Clipping','on',...
'Position',[0.0311750599520384 0.0221565731166913 0.928057553956834 0.556868537666174],...
'CreateFcn', {} );
obj.Handles.RefScanBox = uicontrol(...
'Parent',obj.Handles.ActionPanel,...
'Units','normalized',...
'Callback',{},...
'FontSize',10,...
'Position',[0.036 0.902234636871509 0.399477806788512 0.0754189944134078],...
'String','',...
'Style','checkbox',...
'Tag','RefScanBox',...
'Value',0,...
'Enable','off',...
'CreateFcn', {} );
obj.Handles.RefScanButton = uicontrol(...
'Parent',obj.Handles.ActionPanel,...
'Units','normalized',...
'Callback',@RefScanButtonCallback,...
'FontAngle','oblique',...
'FontWeight','bold',...
'ForegroundColor',[0 0 0],...
'Position',[0.1 0.902234636871509 0.399477806788512 0.0754189944134078],...
'String','Load AM',...
'Tag','RefScanButton',...
'CreateFcn', {} );
 setappdata(obj.Handles.RefScanButton,'Application',obj);
obj.Handles.ExtractTextureBox = uicontrol(...
'Parent',obj.Handles.ActionPanel,...
'Units','normalized',...
'Callback',{},...
'FontSize',10,...
'Position',[0.036 0.82122905027933 0.399477806788512 0.0754189944134078],...
'String','Extract Texture',...
'Style','checkbox',...
'Tag','ExtractTextureBox',...
'CreateFcn', {} );
obj.Handles.StartButton = uicontrol(...
'Parent',obj.Handles.ActionPanel,...
'Units','normalized',...
'Callback',@StartButtonCallback,...
'FontAngle','oblique',...
'FontSize',16,...
'FontWeight','bold',...
'ForegroundColor',[0.392156862745098 0.474509803921569 0.635294117647059],...
'Position',[0.0391644908616188 0.0251396648044693 0.919060052219321 0.201117318435754],...
'String','Start Batch Process',...
'Tag','pushbutton3',...
'CreateFcn', {} );
setappdata(obj.Handles.StartButton,'Application',obj);
end
%% CALLBACKS
function PanelCloseCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         delete(obj);
         close all;
end
function InputPathButtonCallback(hObject,eventdata)
         path = uigetdir(pwd,'Input Folder');
         if path==0, return; end
         cd(path);
         obj = getappdata(hObject,'Application');
         set(obj.Handles.InputPathEdit,'string',path);
end
function OutputPathButtonCallback(hObject,eventdata)
         path = uigetdir(pwd,'Output Folder');
         if path==0, return; end
         cd(path);
         obj = getappdata(hObject,'Application');
         set(obj.Handles.OutputPathEdit,'string',path);
end
function RefScanButtonCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         obj.B.RefScan = meshObj.import;
         if isempty(obj.B.RefScan), return; end
         viewer(obj.B.RefScan);
         set(obj.Handles.RefScanBox,'Value',1);     
end
function StartButtonCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         try
            obj.B.InputFolder = get(obj.Handles.InputPathEdit,'string');
         catch %#ok<CTCH>
             msgbox('Invalid Input Folder','error');
             return;
         end
         obj.B.IncludeSubFolders = get(obj.Handles.IncludeSubFoldersBox,'value');
         obj.B.InputFormat = obj.B.ImportFormats{get(obj.Handles.InputFormatMenu,'value')};
         try 
            obj.B.OutputFolder = get(obj.Handles.OutputPathEdit,'string');
         catch %#ok<CTCH>
             msgbox('Invalid OutputFolder Folder','error');
             return;
         end
         obj.B.Overwrite = get(obj.Handles.OverwriteBox,'value');
         obj.B.OutputFormat = obj.B.ExportFormats{get(obj.Handles.OutputFormatMenu,'value')};
         if isempty(obj.B.RefScan), msgbox('Load Anthropometric Mask first','error'); return; end
         obj.B.ExtractTexture = get(obj.Handles.ExtractTextureBox,'value');
         process(obj.B);
end