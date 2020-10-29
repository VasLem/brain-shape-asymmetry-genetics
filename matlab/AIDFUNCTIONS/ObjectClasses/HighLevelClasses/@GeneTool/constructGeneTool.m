function constructGeneTool(obj)
% GENERATING THE MAIN FIGURE PANEL
  obj.Handles.Panel = figure('Tag', 'Panel', ...
                             'Units', 'characters', ...
                             'Position', [1 5 150 50], ...
                             'Name', 'Janus Morpheus: @ 2013 copyright P. Claes & M. Shriver', ...
                             'MenuBar', 'none', ...
                             'NumberTitle', 'off', ...
                             'Color', get(0,'DefaultUicontrolBackgroundColor'), ...
                             'Resize', 'off',...
                             'CloseRequestFcn',@PanelCloseCallback);
   setappdata(obj.Handles.Panel,'Application',obj);                      
% GENERATING THE TAB PANELS                            
 obj.Handles.TabPanels = uitabpanel('Parent',obj.Handles.Panel,...
                                    'TabPosition','lefttop',...
                                    'Units','normalized',...
                                    'Position',[0,0.10,0.5,0.9],...
                                    'Margins',{[0,-1,1,0],'pixels'},...
                                    'Title',{'SEX & Ancestry', 'SNPs', 'Shape Space'},...
                                    'TitleBackgroundColor',get(0,'DefaultUicontrolBackgroundColor'),...
                                    'TitleForegroundColor',[0 0 0.502],...
                                    'PanelBackgroundColor',get(0,'DefaultUicontrolBackgroundColor'),...
                                    'FrameBackgroundColor',get(0,'DefaultUicontrolBackgroundColor'),...
                                    'PanelBorderType','none',...
                                    'ResizeFcn',{});
tmp = getappdata(obj.Handles.TabPanels,'panels');
obj.Handles.PC.Panel = tmp(3);
constructPanelPC(obj);
obj.Handles.C.Panel = tmp(1);
constructPanelC(obj);
obj.Handles.X.Panel = tmp(2);
constructPanelX(obj);
% % GENERATING THE BUTTONS
obj.Handles.RenderMode = uicontrol(...
                                'Parent',obj.Handles.Panel,...
                                'Units','normalized',...
                                'BackgroundColor',[1 1 1],...
                                'Callback',@RenderModeMenuCallback,...
                                'Position',[0.65 0.10 0.15 0.05],...
                                'String',{'shape' 'appearance' 'curvature' 'area' 'normal'},...
                                'Style','popupmenu',...
                                'Value',1,...
                                'TooltipString','Select RenderMode',...
                                'CreateFcn',{},...
                                'Tag','RenderMode Menu');
setappdata(obj.Handles.RenderMode,'Application',obj);
obj.Handles.UpdateButton = uicontrol('Parent', obj.Handles.Panel, ...
                                   'Tag', 'Update Button', ...
                                   'UserData', [], ...
                                   'Style', 'pushbutton', ...
                                   'Units', 'normalized', ...
                                   'Position', [0.01 0.02 0.15 0.05], ...
                                   'FontWeight', 'bold', ...
                                   'ForegroundColor', [0 0 0.502], ...
                                   'String', 'Update FChar',...
                                   'Callback', @UpdateButtonCallback);
setappdata(obj.Handles.UpdateButton,'Application',obj);
obj.Handles.AvgButton = uicontrol('Parent', obj.Handles.Panel, ...
                                   'Tag', 'Avg Button', ...
                                   'UserData', [], ...
                                   'Style', 'pushbutton', ...
                                   'Units', 'normalized', ...
                                   'Position', [0.17 0.02 0.15 0.05], ...
                                   'FontWeight', 'bold', ...
                                   'ForegroundColor', [0 0 0.502], ...
                                   'String', 'Average',...
                                   'Callback', @AvgButtonCallback);
setappdata(obj.Handles.AvgButton,'Application',obj);
obj.Handles.SetButton = uicontrol('Parent', obj.Handles.Panel, ...
                                   'Tag', 'Set Button', ...
                                   'UserData', [], ...
                                   'Style', 'pushbutton', ...
                                   'Units', 'normalized', ...
                                   'Position', [0.33 0.02 0.15 0.05], ...
                                   'FontWeight', 'bold', ...
                                   'ForegroundColor', [0 0 0.502], ...
                                   'String', 'Set',...
                                   'Callback', @SetButtonCallback);
setappdata(obj.Handles.SetButton,'Application',obj);
obj.Handles.ResetButton = uicontrol('Parent', obj.Handles.Panel, ...
                                   'Tag', 'Reset Button', ...
                                   'UserData', [], ...
                                   'Style', 'pushbutton', ...
                                   'Units', 'normalized', ...
                                   'Position', [0.49 0.02 0.15 0.05], ...
                                   'FontWeight', 'bold', ...
                                   'ForegroundColor', [0 0 0.502], ...
                                   'String', 'Reset',...
                                   'Callback', @ResetButtonCallback);
setappdata(obj.Handles.ResetButton,'Application',obj);
obj.Handles.ImportButton = uicontrol('Parent', obj.Handles.Panel, ...
                                    'Tag', 'Import Button', ...
                                	'UserData', [], ...
                                	'Style', 'pushbutton', ...
                                    'Units', 'normalized', ...
                                	'Position', [0.65 0.02 0.15 0.05], ...
                                	'FontWeight', 'bold', ...
                                	'ForegroundColor', [0 0 0.502], ...
                                	'String', 'Import Face',...
                                    'Enable','on',...
                                    'Callback', @ImportButtonCallback);
setappdata(obj.Handles.ImportButton,'Application',obj);
% obj.Handles.LoadModelButton = uicontrol('Parent', obj.Handles.Panel, ...
%                                     'Tag', 'Load Model Button', ...
%                                 	'UserData', [], ...
%                                 	'Style', 'pushbutton', ...
%                                     'Units', 'normalized', ...
%                                 	'Position', [0.81 0.02 0.15 0.05], ...
%                                 	'FontWeight', 'bold', ...
%                                 	'ForegroundColor', [0 0 0.502], ...
%                                 	'String', 'Load Model',...
%                                     'Enable','off',...
%                                     'Callback', @LoadModelButtonCallback);
% setappdata(obj.Handles.LoadModelButton,'Application',obj);
obj.Handles.ExportButton = uicontrol('Parent', obj.Handles.Panel, ...
                                    'Tag', 'Export Button', ...
                                	'UserData', [], ...
                                	'Style', 'pushbutton', ...
                                    'Units', 'normalized', ...
                                	'Position', [0.81 0.02 0.15 0.05], ...
                                	'FontWeight', 'bold', ...
                                	'ForegroundColor', [0 0 0.502], ...
                                	'String', 'Export Face',...
                                    'Enable','on',...
                                    'Callback', @ExportButtonCallback);
setappdata(obj.Handles.ExportButton,'Application',obj);
obj.Handles.PC.Axes =  axes('Parent', obj.Handles.Panel, ...
                                'Tag', 'render_axes', ...
                                'UserData', obj, ...
                                'Units', 'normalized', ...
                                'Position',[0.55 0.20 0.4 0.7],...
                                'XColor',[0 0 0],...
                                'YColor',[0 0 0],...
                                'ZColor',[0 0 0]);
        view(obj.Handles.PC.Axes,-270,90);
        setappdata(obj.Handles.PC.Axes,'Application',obj);
end
%% MAIN PANEL CALLBACKS
function PanelCloseCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         delete(obj);
         close all;
end
function ResetButtonCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         updateCurrentShape(obj,'RESET');
end
function SetButtonCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         updateCurrentShape(obj,'SET');
end
function ImportButtonCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         updateCurrentShape(obj,'IMPORT');
end
function ExportButtonCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         [filename, path, index] = uiputfile({'*.xls';'*.xlsx'}, 'Face Export To Excell');
         cd(path)
         filename = filename(1:end-4);
         headings = cell(1,3);
         headings{1} = 'PC';
         headings{2} = 'C';
         headings{3} = 'X';
         xlswrite(filename,headings,1,'A1');
         xlswrite(filename,obj.CurrentCoeff',1,'A2');
         xlswrite(filename,obj.CurrentC',1,'B2');
         xlswrite(filename,obj.CurrentX',1,'C2');
end
function UpdateButtonCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         updateFC(obj);
end
function AvgButtonCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         updateCurrentShape(obj,'AVERAGE');
end
function LoadModelButtonCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         out = getMatFile('select shape model');
         if isempty(out), return; end
         switch out{1}.Type
             case 'CombinedPLSRConditionedShapeModel'
                 obj.PredModel = CombinedPLSRConditionedShapeModel(out{1});
             otherwise
                 error('Wrong model input');
         end
end
function RenderModeMenuCallback(hObject,eventdata)
          obj = getappdata(hObject,'Application');
          list = get(obj.Handles.RenderMode,'String');
          obj.RenderMode = list{get(obj.Handles.RenderMode,'Value')};
end
%% PC Panel creation functions
function constructPanelPC(obj)
          obj.Handles.PC.Table = uitable(...
                                            'Parent',obj.Handles.PC.Panel,...
                                            'Units','normalized',...
                                            'BackgroundColor',[1 1 1;0.96078431372549 0.96078431372549 0.96078431372549],...
                                            'ColumnName',{'Sel.', 'PC', 'Value', 'mean', 'std'},...
                                            'ColumnFormat',{'logical','char','numeric','numeric','numeric'},...
                                            'ColumnEditable',[true false true false false],...
                                            'ColumnWidth',{40 40 60 60 60},...
                                            'Data',{},...
                                            'Position',[0.01 0.11 0.80 0.85],...
                                            'UserData',[],...
                                            'Tag','PCTable',...
                                            'CellEditCallback', @PCTableCallback,...
                                            'RowName', []);
          setappdata(obj.Handles.PC.Table,'Application',obj);                                        
          obj.Handles.PC.Slider = uicontrol(...
                                            'Parent',obj.Handles.PC.Panel,...
                                            'Units','normalized',...
                                            'BackgroundColor',[0.9 0.9 0.9],...
                                            'Callback',@PCSliderCallback,...
                                            'Position',[0.01 0.01 0.98 0.08],...
                                            'String',{  'Slider' },...
                                            'Style','slider',...
                                            'CreateFcn',{},...
                                            'Max',obj.PCScale,...
                                            'Min',-obj.PCScale,...
                                            'Tag','PC Slider');
         setappdata(obj.Handles.PC.Slider,'Application',obj);
         obj.Handles.PC.Scale = uicontrol(...
                                'Parent',obj.Handles.PC.Panel,...
                                'Units','normalized',...
                                'BackgroundColor',[1 1 1],...
                                'Callback',@PCScaleMenuCallback,...
                                'Position',[0.85 0.05 0.15 0.1],...
                                'String',{'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'},...
                                'Style','popupmenu',...
                                'Value',obj.PCScale,...
                                'TooltipString','Select Slider Scale',...
                                'CreateFcn',{},...
                                'Tag','PC Scale Menu');
        setappdata(obj.Handles.PC.Scale,'Application',obj);
end
%% PC PANEL CALLBACKS
function PCSliderCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         value = get(obj.Handles.PC.Slider,'Value');
         maximum = get(obj.Handles.PC.Slider,'max');
         minimum = get(obj.Handles.PC.Slider,'min');
         if value > maximum
            value = maximum;
            set(obj.Handles.PC.Slider,'Value',value);
         end
         if value < minimum;
            value = minimum;
            set(obj.Handles.PC.Slider,'Value',value);
         end
         updateShapeTable(obj,obj.ActivePC,value);
         updateCurrentShape(obj,'PC');
end
function PCScaleMenuCallback(hObject,eventdata)
          obj = getappdata(hObject,'Application');
          list = get(obj.Handles.PC.Scale,'String');
          obj.PCScale = str2double(list{get(obj.Handles.PC.Scale,'Value')});
end
function PCTableCallback(hObject,eventdata)
 obj = getappdata(hObject,'Application');
 if isempty(obj.ShapeModel), return; end
 row = eventdata.Indices(1);
 Column = eventdata.Indices(2);
 new = eventdata.NewData;
 switch Column
     case 1% change of active scan
         if obj.ActivePC == row  
             obj.ActivePC = 0;
         else
             if ~(obj.ActivePC==0)
                 Data = get(obj.Handles.PC.Table,'Data');
                 Data{obj.ActivePC,1} = false;
                 set(obj.Handles.PC.Table,'Data',Data);
             end
             obj.ActivePC = row;
         end
     case 3% change of PC value
         Data = get(obj.Handles.PC.Table,'Data');
         maximum = Data{row,4}+obj.PCScale*Data{row,5};
         minimum = Data{row,4}-obj.PCScale*Data{row,5};
         if new > maximum
            new = maximum;   
         end
         if new < minimum;
            new = minimum;
         end
         Data{row,Column} = new;
         if ~(obj.ActivePC==row)          
             if ~(obj.ActivePC == 0)
                Data{obj.ActivePC,1} = false;
             end
             Data{row,1} = true;
         end
         set(obj.Handles.PC.Table,'Data',Data);
         updateCurrentShape(obj,'PC');
         obj.ActivePC = row;
     otherwise
         % do nothing
 end
end
%% C Panel creation functions
function constructPanelC(obj)
          obj.Handles.C.Table = uitable(...
                                            'Parent',obj.Handles.C.Panel,...
                                            'Units','normalized',...
                                            'BackgroundColor',[1 1 1;0.96078431372549 0.96078431372549 0.96078431372549],...
                                            'ColumnName',{'Sel.', 'Pred.', 'Value', 'mean', 'std'},...
                                            'ColumnFormat',{'logical','char','numeric','numeric','numeric'},...
                                            'ColumnEditable',[true false true false false],...
                                            'ColumnWidth',{40 40 60 60 60},...
                                            'Data',{},...
                                            'Position',[0.01 0.11 0.80 0.85],...
                                            'UserData',[],...
                                            'Tag','PCTable',...
                                            'CellEditCallback', @CTableCallback,...
                                            'RowName', []);
          setappdata(obj.Handles.C.Table,'Application',obj);                                        
          obj.Handles.C.Slider = uicontrol(...
                                            'Parent',obj.Handles.C.Panel,...
                                            'Units','normalized',...
                                            'BackgroundColor',[0.9 0.9 0.9],...
                                            'Callback',@CSliderCallback,...
                                            'Position',[0.01 0.01 0.98 0.08],...
                                            'String',{  'Slider' },...
                                            'Style','slider',...
                                            'CreateFcn',{},...
                                            'Max',obj.CScale,...
                                            'Min',-obj.CScale,...
                                            'Tag','C Slider');
         setappdata(obj.Handles.C.Slider,'Application',obj);
         obj.Handles.C.Scale = uicontrol(...
                                'Parent',obj.Handles.C.Panel,...
                                'Units','normalized',...
                                'BackgroundColor',[1 1 1],...
                                'Callback',@CScaleMenuCallback,...
                                'Position',[0.85 0.05 0.15 0.1],...
                                'String',{'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'},...
                                'Style','popupmenu',...
                                'Value',obj.CScale,...
                                'TooltipString','Select Slider Scale',...
                                'CreateFcn',{},...
                                'Tag','C Scale Menu');
        setappdata(obj.Handles.C.Scale,'Application',obj);
end
%% C PANEL CALLBACKS
function CSliderCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         value = get(obj.Handles.C.Slider,'Value');
         maximum = get(obj.Handles.C.Slider,'max');
         minimum = get(obj.Handles.C.Slider,'min');
         if value > maximum
            value = maximum;
            set(obj.Handles.C.Slider,'Value',value);
         end
         if value < minimum;
            value = minimum;
            set(obj.Handles.C.Slider,'Value',value);
         end
         updateCTable(obj,obj.ActiveC,value);
         updateCurrentShape(obj,'C');
end
function CScaleMenuCallback(hObject,eventdata)
          obj = getappdata(hObject,'Application');
          list = get(obj.Handles.C.Scale,'String');
          obj.CScale = str2double(list{get(obj.Handles.C.Scale,'Value')});
end
function CTableCallback(hObject,eventdata)
 obj = getappdata(hObject,'Application');
 if isempty(obj.PredModel), return; end
 row = eventdata.Indices(1);
 Column = eventdata.Indices(2);
 new = eventdata.NewData;
 switch Column
     case 1% change of active scan
         if obj.ActiveC == row  
             obj.ActiveC = 0;
         else
             if ~(obj.ActiveC==0)
                 Data = get(obj.Handles.C.Table,'Data');
                 Data{obj.ActiveC,1} = false;
                 set(obj.Handles.C.Table,'Data',Data);
             end
             obj.ActiveC = row;
         end
     case 3 % change of C value
         Data = get(obj.Handles.C.Table,'Data');
         maximum = Data{row,4}+obj.CScale*Data{row,5};
         minimum = Data{row,4}-obj.CScale*Data{row,5};
         if new > maximum
            new = maximum;   
         end
         if new < minimum;
            new = minimum;
         end
         Data{row,Column} = new;
         if ~(obj.ActiveC==row)          
             if ~(obj.ActiveC == 0)
                Data{obj.ActiveC,1} = false;
             end
             Data{row,1} = true;
         end
         set(obj.Handles.C.Table,'Data',Data);
         updateCurrentShape(obj,'C');
         obj.ActiveC = row;
     otherwise
         % do nothing
 end
end
%% X Panel creation functions
function constructPanelX(obj)
          obj.Handles.X.Table = uitable(...
                                            'Parent',obj.Handles.X.Panel,...
                                            'Units','normalized',...
                                            'BackgroundColor',[1 1 1;0.96078431372549 0.96078431372549 0.96078431372549],...
                                            'ColumnName',{'Sel.', 'Pred.', 'Value', 'mean', 'std'},...
                                            'ColumnFormat',{'logical','char','numeric','numeric','numeric'},...
                                            'ColumnEditable',[true false true false false],...
                                            'ColumnWidth',{40 40 60 60 60},...
                                            'Data',{},...
                                            'Position',[0.01 0.11 0.80 0.85],...
                                            'UserData',[],...
                                            'Tag','PCTable',...
                                            'CellEditCallback', @XTableCallback,...
                                            'RowName', []);
          setappdata(obj.Handles.X.Table,'Application',obj);                                        
          obj.Handles.X.Slider = uicontrol(...
                                            'Parent',obj.Handles.X.Panel,...
                                            'Units','normalized',...
                                            'BackgroundColor',[0.9 0.9 0.9],...
                                            'Callback',@XSliderCallback,...
                                            'Position',[0.01 0.01 0.98 0.08],...
                                            'String',{  'Slider' },...
                                            'Style','slider',...
                                            'CreateFcn',{},...
                                            'Max',obj.XScale,...
                                            'Min',-obj.XScale,...
                                            'Tag','X Slider');
         setappdata(obj.Handles.X.Slider,'Application',obj);
         obj.Handles.X.Scale = uicontrol(...
                                'Parent',obj.Handles.X.Panel,...
                                'Units','normalized',...
                                'BackgroundColor',[1 1 1],...
                                'Callback',@XScaleMenuCallback,...
                                'Position',[0.85 0.05 0.15 0.1],...
                                'String',{'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'},...
                                'Style','popupmenu',...
                                'Value',obj.XScale,...
                                'TooltipString','Select Slider Scale',...
                                'CreateFcn',{},...
                                'Tag','X Scale Menu');
        setappdata(obj.Handles.X.Scale,'Application',obj);
end
%% PC PANEL CALLBACKS
function XSliderCallback(hObject,eventdata)
         obj = getappdata(hObject,'Application');
         value = get(obj.Handles.X.Slider,'Value');
         maximum = get(obj.Handles.X.Slider,'max');
         minimum = get(obj.Handles.X.Slider,'min');
         if value > maximum
            value = maximum;
            set(obj.Handles.X.Slider,'Value',value);
         end
         if value < minimum;
            value = minimum;
            set(obj.Handles.X.Slider,'Value',value);
         end
         updateXTable(obj,obj.ActiveX,value);
         updateCurrentShape(obj,'X');
end
function XScaleMenuCallback(hObject,eventdata)
          obj = getappdata(hObject,'Application');
          list = get(obj.Handles.X.Scale,'String');
          obj.XScale = str2double(list{get(obj.Handles.X.Scale,'Value')});
end
function XTableCallback(hObject,eventdata)
 obj = getappdata(hObject,'Application');
 if isempty(obj.PredModel), return; end
 row = eventdata.Indices(1);
 Column = eventdata.Indices(2);
 new = eventdata.NewData;
 switch Column
     case 1% change of active scan
         if obj.ActiveX == row  
             obj.ActiveX = 0;
         else
             if ~(obj.ActiveX==0)
                 Data = get(obj.Handles.X.Table,'Data');
                 Data{obj.ActiveX,1} = false;
                 set(obj.Handles.X.Table,'Data',Data);
             end
             obj.ActiveX = row;
         end
     case 3
         Data = get(obj.Handles.X.Table,'Data');
         maximum = Data{row,4}+obj.XScale*Data{row,5};
         minimum = Data{row,4}-obj.XScale*Data{row,5};
         if new > maximum
            new = maximum;   
         end
         if new < minimum;
            new = minimum;
         end
         Data{row,Column} = new;
         if ~(obj.ActiveX==row)          
             if ~(obj.ActiveX == 0)
                Data{obj.ActiveX,1} = false;
             end
             Data{row,1} = true;
         end
         set(obj.Handles.X.Table,'Data',Data);
         obj.ActiveX = row;
         updateCurrentShape(obj,'X');
     otherwise
         % do nothing
 end
end
