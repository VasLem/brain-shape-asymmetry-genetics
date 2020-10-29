function constructViewer(obj)
%% Creation
% --- FIGURE CREATION -------------------------------------
obj.Figure = figure(	'Tag', '3D Viewer', ...
	'Units', 'characters', ...
	'Position', [103 24 100 30], ...
	'Name', obj.Tag, ...
	'MenuBar', 'none', ...
    'renderer','zbuffer',...
	'NumberTitle', 'off', ...
	'Color', [0    0    0], ...
	'Resize', 'on',...
    'CloseRequestFcn',@figure_close_Callback,...
    'WindowButtonDownFcn',@BtnDown,...
    'WindowButtonMotionFcn', @BtnMotion,...
    'WindowButtonUpFcn',@BtnUp,...
    'WindowKeyPressFcn',@KeyPress,...
    'WindowKeyReleaseFcn',@KeyRelease,...
    'WindowScrollWheelFcn',@ScrollWeel,...
    'UserData',obj);
% warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
% jframe=get(obj.Figure,'javaframe');
% tmp = which('constructViewer');
% tmp = [tmp(1:end-17) 'unknownicon.gif'];
% jIcon=javax.swing.ImageIcon(tmp);
% jframe.setFigureIcon(jIcon);

% --- AXES CREATION -------------------------------------
obj.RenderAxes = axes(	'Parent', obj.Figure, ...
	'Tag', 'render_axes', ...
	'UserData', obj, ...
	'Units', 'normalized', ...
    'Position',[0.0946428571428572 0.0761904761904762 0.805357142857143 0.835714285714286],...
    'Color',[0 0 0],...
    'XColor',[1 1 1],...
    'YColor',[1 1 1],...
    'ZColor',[1 1 1]);
xlabel(obj.RenderAxes,'X');
ylabel(obj.RenderAxes,'Y');
zlabel(obj.RenderAxes,'Z');
axis(obj.RenderAxes,'off','equal');
hold(obj.RenderAxes,'on');
% camproj(obj.RenderAxes,'orthographic');
% --- SCENE LIGHT CONTRUCTION ------------------------------
obj.SceneLight = camlight('headlight','infinite');
set(obj.SceneLight,'Visible','off');
% --- TOOLBAR CREATION -------------------------------------
constructViewerToolbar(obj);
% --- CONTEXT MENU CREATION --------------------------------
constructContextMenu(obj);
set(obj.Figure,'UIContextMenu',obj.ContextMenu.h);
% set(obj.RenderAxes,'UIContextMenu',obj.ContextMenu.h);
obj.CamProjection = 'orthographic';
% --- FILTER SETTING GUI CREATION ----
constructFilterSettingsWindow(obj);


end

%% CALLBACKS

function figure_close_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         obj.Visible = false;
         if isempty(obj.Parent)
             delete(obj);
         else
             handles = guidata(obj.Parent);
             handles.viewerChange(handles,'close');
         end
end

function KeyPress(varargin)
 obj = get(varargin{1},'UserData');
 if strcmp(obj.Status,'Busy'), return; end
 %varargin{2}
 switch obj.SelectionMode
     case 'landmark'
         switch varargin{2}.Key
             case 'space'
                  if isempty(obj.CurrentMesh), 
                     b = msgbox(obj.NoCurrentMeshMsg,'Warning');
                     waitfor(b);
                     return; 
                  end
                  screenPoint2currentPoint(obj)
                  obj.ActiveKey = 'space';
                  set(obj.Figure,'WindowKeyPressFcn',@dummyKeyPress);
                  set(obj.Figure,'UIContextMenu',[]);                 
             case 'backspace'
                 updateLMSelection(obj,'Clear'); 
                 obj.ActiveKey = 'backspace';
             case 'p'
                 updateLMSelection(obj,'Link Pose LM');
                 obj.ActiveKey = 'p';
                 if strcmp(obj.CalledProcess,'indicatePoseLM')
                         delete(obj);return;
                 end
             case 'c'
                 updateLMSelection(obj,'Link Custom LM');
                 obj.ActiveKey = 'c';
             case 'u'
                 try
                   obj.UserData = obj.LandmarkSelection.Vertices;
                   obj.UserData
                   obj.ActiveKey = 'u';
                   obj.CurrentMesh.UserData = clone(obj.LandmarkSelection);
                 catch
                 end
             otherwise
                 otherKeys(varargin{2}.Key,obj);
         end
     case {'area' 'brush' 'fill'}
         switch varargin{2}.Key
             case 'space'
                 if isempty(obj.CurrentMesh), uiwait(msgbox(obj.NoCurrentMeshMsg,'Warning','modal')); return; end
                 obj.ActiveKey = 'space';
                 linkAreaSelection(obj);
                 set(obj.Figure,'WindowKeyPressFcn',@dummyKeyPress);
                 set(obj.Figure,'UIContextMenu',[]);
                 % special for brush mode
                 if strcmp(obj.SelectionMode,'brush')
                     linkAreaSelectionRadius(obj); 
                     obj.Action = 'Brush Dummy';
                     updateAreaSelectionRadius(obj,'Add Area', select3DRadiusArea(obj));
                     set(obj.Figure,'WindowButtonMotionFcn',@BtnMotion);
                 end
             case 'delete'
                 tmp = varargin{2}.Modifier;
                 if isempty(tmp)
                    mode = 'normal';
                 else
                    mode = tmp{1};
                 end
                 switch mode
                     case 'alt'
                         updateAreaSelection(obj,'Crop');
                     otherwise
                         updateAreaSelection(obj,'Delete');
                 end
             case 'backspace'
                 updateAreaSelection(obj,'Clear');
             case 'i'
                 updateAreaSelection(obj,'Invert');
             case 's'
                 updateAreaSelection(obj,'Smooth Surface','Busy');
             case 'c'
                 updateAreaSelection(obj,'Smooth Color','Busy');
             case {'add' 'equal'}
                 updateAreaSelection(obj,'Subdivide Triangles');
             case {'subtract' 'hyphen'}
                 updateAreaSelection(obj,'Reduce Triangles');
             case 'u'
                 try
                    obj.UserData = obj.AreaSelection.VerticesIndex;
                    obj.UserData'
                    obj.ActiveKey = 'u';
                    obj.CurrentMesh.UserData = obj.UserData;
                 catch 
                 end
             otherwise
                 otherKeys(varargin{2}.Key,obj);
         end
     case 'full'
         switch varargin{2}.Key
             case 'space'
                 obj.ActiveKey = 'space';
                 set(obj.Figure,'WindowKeyPressFcn',@dummyKeyPress);
                 set(obj.Figure,'UIContextMenu',[]);
             case 'delete'
                 tmp = varargin{2}.Modifier;
                 if isempty(tmp)
                    mode = 'normal';
                 else
                    mode = tmp{1};
                 end
                 switch mode
                     case 'alt'
                         updateFullScanSelection(obj,'Crop');
                     otherwise
                         updateFullScanSelection(obj,'Delete');
                 end
             case 'backspace'
                 updateFullScanSelection(obj,'Clear');
             case 'i'
                 updateFullScanSelection(obj,'Invert');
             case 's'
                 updateFullScanSelection(obj,'Smooth Surface');
             case 'c'
                 updateFullScanSelection(obj,'Smooth Color');
             case {'add' 'equal'}
                 updateFullScanSelection(obj,'Subdivide Triangles');
             case {'subtract' 'hyphen'}
                 updateFullScanSelection(obj,'Reduce Triangles','Busy');
             otherwise
                 otherKeys(varargin{2}.Key,obj);
         end
     otherwise
         otherKeys(varargin{2}.Key,obj);
 end
 try
    if obj.Record && ~strcmp(obj.Status,'Recording')
       captureFrame(obj);
    end
 catch
 end
    
end

function otherKeys(key,obj)
% key
       switch key
           case 'l'
               switch obj.SceneLightVisible
                   case true
                      set(obj.Toolbar.light_toggle,'State','off');
                   case false
                      set(obj.Toolbar.light_toggle,'State','on'); 
               end
           case 'f1'
               updateFullScanSelection(obj,'Rendering Surface','Solid');
           case 'f2'
               updateFullScanSelection(obj,'Rendering Surface','Wireframe');
           case 'f3'
               updateFullScanSelection(obj,'Rendering Surface','Points');
           case 'f4'
               updateFullScanSelection(obj,'Rendering Surface','Solid/Wireframe');
           case 'f5'
               updateFullScanSelection(obj,'Rendering Color','Texture');
           case 'f6'
               updateFullScanSelection(obj,'Rendering Color','Indexed');
           case 'f7'
               updateFullScanSelection(obj,'Rendering Color','Single');
           case 'f8'
               updateFullScanSelection(obj,'Rendering Material','Facial');
           case 'f9'
               updateFullScanSelection(obj,'Rendering Material','Dull');
           case 'f10'
               updateFullScanSelection(obj,'Rendering Material','Shiny');
           case 'f11'
               updateFullScanSelection(obj,'Rendering Material','Metal');
           case 'f12'
               updateFullScanSelection(obj,'Rendering Material','Default');
           case {'0' '1' '2' '3' '4' '5' '6' '7' '8' '9' 'numpad0' 'numpad1' 'numpad2' 'numpad3' 'numpad4' 'numpad5' 'numpad6' 'numpad7' 'numpad8' 'numpad9'}
                updateFullScanSelection(obj,'Rendering Transparancy',(10-str2double(key(end)))/10);
           case 'return'
               if strcmp(obj.CalledProcess,'Batch Processing')
                         delete(obj);return;
               end
           otherwise
               return
       end
end

function dummyKeyPress(varargin) 
     % Do nothing
end

function KeyRelease(varargin)
         obj = get(varargin{1},'UserData');        
         switch obj.SelectionMode
             case {'landmark' 'area' 'fill' 'brush' 'full'}
                 switch varargin{2}.Key
                     case 'space'
                         set(obj.Figure,'WindowKeyPressFcn',@KeyPress);
                         set(obj.Figure,'UIContextMenu',obj.ContextMenu.h);
                         % special for brush mode
                         if strcmp(obj.SelectionMode,'brush') 
                             updateAreaSelectionRadius(obj,'Clear');
                             obj.Action = 'none';
                             set(obj.Figure,'WindowButtonMotionFcn',[]);
                         end
                     otherwise
                 end
             otherwise
         end
         obj.ActiveKey = 'none';
         %disp('Key release')
         if obj.Record
            captureFrame(obj);
         end
end

function BtnDown(varargin)
     obj = get(varargin{1},'UserData');
     pt = hgconvertunits(obj.Figure,[0 0 get(obj.Figure,'CurrentPoint')],get(obj.Figure,'Units'),'pixels',0);
     obj.MotionData.CurrentPoint = pt(3:4);
     switch obj.ActiveKey
         case 'space'
             BtnDownSelectionMode(obj);
         otherwise
             if strcmp(obj.Mode,'none'), return; end
             BtnDownMode(obj);
     end
     if obj.Record && ~strcmp(obj.Status,'Recording')
       captureFrame(obj);
     end
end

function BtnDownSelectionMode(obj)
     switch obj.SelectionMode
        case 'landmark'           
             [p, v, vi, f, fi]  = select3DPoint(obj);
             if isempty(p), return; end
             switch get(obj.Figure, 'selectiontype')
                 case 'normal'
                      updateLMSelection(obj,'Add LM',p,vi,fi);
                 case 'alt'
                      updateLMSelection(obj,'Delete LM',p);
             end
        case 'area'
             switch get(obj.Figure, 'selectiontype')
                 case 'normal'% Box add selection
                      linkAreaSelection(obj);
                      areaindex = select3DRbBoxArea(obj);
                      if ~isempty(areaindex), updateAreaSelection(obj,'Add Area',areaindex); end
                 case 'alt'% Box remove selection
                      linkAreaSelection(obj);
                      areaindex = select3DRbBoxArea(obj);
                      if ~isempty(areaindex), updateAreaSelection(obj,'Remove Area',areaindex); end
                 case 'extend'
                      updateAreaSelection(obj,'Invert');
                 case 'open'
                      updateAreaSelection(obj,'Clear');
                 otherwise
                     return;
             end         
        case 'fill'
            switch get(obj.Figure, 'selectiontype')
                 case 'normal'% Connected add selection
                      obj.Status = 'Busy';drawnow;
                      linkAreaSelection(obj);
                      areaindex = select3DConnectedArea(obj);
                      if ~isempty(areaindex), updateAreaSelection(obj,'Add Area',areaindex); end
                      obj.Status = 'Ready';
                 case 'alt'% Connected remove selection
                      obj.Status = 'Busy';drawnow;
                      linkAreaSelection(obj);
                      areaindex = select3DConnectedArea(obj);
                      if ~isempty(areaindex), updateAreaSelection(obj,'Remove Area',areaindex); end
                      obj.Status = 'Ready';
                 case 'extend'
                     updateAreaSelection(obj,'Invert'); 
                 case 'open'
                     updateAreaSelection(obj,'Clear');
                 otherwise
                     return;
            end      
        case 'brush'
             switch get(obj.Figure, 'selectiontype')
                 case 'normal'% Radius add selection
                      updateAreaSelectionRadius(obj,'Clear');
                      linkAreaSelection(obj);
                      areaindex = select3DRadiusArea(obj);
                      if ~isempty(areaindex), updateAreaSelection(obj,'Add Area',areaindex); end
                      obj.Action = 'Brush Selection';
                      set(obj.Figure,'WindowButtonMotionFcn',@BtnMotion);
                 case 'alt'% Radius remove selection
                      linkAreaSelection(obj);
                      areaindex = select3DRadiusArea(obj);
                      if ~isempty(areaindex), updateAreaSelection(obj,'Remove Area',areaindex); end
                      obj.Action = 'Brush Deselection';
                      obj.AreaSelectionRadius.SingleColor = [0 0 0.502];
                      set(obj.Figure,'WindowButtonMotionFcn',@BtnMotion);                      
                 case 'extend'
                      updateAreaSelection(obj,'Invert');
                 case 'open'
                      updateAreaSelection(obj,'Clear');
                 otherwise
                     return;
             end
         case 'full'
             selection = get(obj.Figure, 'CurrentObject');
             if ~strcmp(get(selection,'Type'),'patch'), return; end
             meshobject = get(selection,'UserData');
             if ~strcmp(class(meshobject),'meshObj'), return; end
             if~validH(meshobject,'Border'), border(meshobject); end
             switch get(obj.Figure, 'Selectiontype');
                 case 'normal'% full scan Selection
                     if meshobject.Selected                    
                         meshobject.Border.SingleColor = [0.8 0 0];
                         pause(0.1);
                     end
                     meshobject.Selected = true;
                     obj.Action = 'Full Scan Selected';
                     if isempty(obj.Parent), return; end
                     % if called with Scan3D feedback to main GUI
                     handles = guidata(obj.Parent);
                     handles.viewerChange(handles,'Scan Selected');                                 
                 case 'alt'% full scan deselection
                     meshobject.Selected = false;
                     if isempty(obj.Parent), return; end
                     % if called with Scan3D feedback to main GUI
                     handles = guidata(obj.Parent);
                     handles.viewerChange(handles,'Scan Selected');           
                 case 'extend'% ?
                 case 'open'% ? 
                 otherwise
             end
        otherwise
            return
    end
end

function BtnDownMode(obj)
      switch obj.Mode
         case 'camera'
             axis(obj.RenderAxes,'vis3d');
             switch get(obj.Figure, 'selectiontype')
                 case 'normal'
                     obj.Action = 'rotate camera';
                 case 'alt'
                     obj.Action = 'pan camera';
                 case 'extend'
                     obj.Action = 'zoom camera';
                 case 'open'
                     resetCamera(obj);
                     return;
                 otherwise
                     return;
             end
         case 'light'
             axis(obj.RenderAxes,'vis3d');
             if ~obj.SceneLightVisible, return; end
             switch get(obj.Figure, 'selectiontype')
                 case 'normal'% free light rotation
                     obj.Action = 'rotate light';
                 case 'alt'% change light mode
                     switch obj.SceneLightMode
                         case 'infinite'
                             obj.SceneLightMode = 'local';
                         case 'local'
                             obj.SceneLightMode = 'infinite';                         
                         otherwise
                             return;
                     end
                     return;
                 case 'extend'
                     switch obj.SceneLightRotMode;
                         case 'vertical'
                             obj.SceneLightRotMode = 'horizontal';
                         case 'horizontal'
                             obj.SceneLightRotMode = 'vertical';
                         otherwise
                             return;
                     end
                     return
                 case 'open'
                     resetLight(obj);
                     return;
                 otherwise
                     return;
             end
         otherwise
             return;
     end
     set(obj.Figure,'WindowButtonMotionFcn',@BtnMotion);
end

function BtnUp(varargin)
         obj = get(varargin{1},'UserData');
         switch obj.Action
             case 'Full Scan Selected'
                 if strcmp(obj.ActiveKey,'space')
                    scrn_pt = get(0, 'PointerLocation');
                    if ~PointerOnFigure(obj.Figure,scrn_pt)
                       figure_list = findobj('Type','figure');
                       for f=1:1:length(figure_list)
                           if ~PointerOnFigure(figure_list(f),scrn_pt), continue; end
                           UD = get(figure_list(f),'UserData');
                           if ~strcmp(class(UD),'viewer3DObj'), continue; end
                           if (UD==obj), break; end
                           selection = get(obj.Figure, 'CurrentObject');
                           meshObject = get(selection,'UserData');
                           meshObject.Axes = UD.RenderAxes;
                           UD.ActiveKey = 'none';
                           break;
                       end
                    end                        
                 end
             otherwise
         end
         % special for active brush mode
         if strcmp(obj.SelectionMode,'brush')&&strcmp(obj.ActiveKey,'space')
            obj.Action = 'Brush Dummy';
            obj.AreaSelectionRadius.SingleColor = [0 0.7 0];
            updateAreaSelectionRadius(obj,'Add Area', select3DRadiusArea(obj));
         else
             obj.Action = 'none';
             set(obj.Figure,'WindowButtonMotionFcn','');
             %set(obj.Figure,'WindowButtonMotionFcn',@dummyBtnMotion);
         end   
         drawnow;
         try
            if obj.Record && ~strcmp(obj.Status,'Recording')
               captureFrame(obj);
            end
         catch
         end
end

function dummyBtnMotion(varargin)% does not work propperly
         obj = get(varargin{1},'UserData');
         if obj.Record
            pause(0.1);
            captureFrame(obj);
         end
end

function BtnMotion(varargin)
         obj = get(varargin{1},'UserData');
         if strcmp(obj.Status,'Recording'), return; end
         if strcmp(obj.Action,'none'), return; end
         pt = hgconvertunits(obj.Figure,[0 0 get(obj.Figure,'CurrentPoint')],get(obj.Figure,'Units'),'pixels',0);
         CurrentPoint = pt(3:4);
         deltaPix  = CurrentPoint-obj.MotionData.CurrentPoint;
         obj.MotionData.CurrentPoint = CurrentPoint;
         switch obj.Action
             case 'rotate camera'
                 rotateCamera(obj,deltaPix);
             case 'pan camera'
                 panCamera(obj,deltaPix);
             case 'zoom camera'
                 zoomCamera(obj,deltaPix,[]);
             case 'rotate light'
                 rotateLight(obj,deltaPix);
             case 'Brush Selection'
                 if ~strcmp(obj.ActiveKey,'space'), obj.Action = 'none'; return; end
                 updateAreaSelection(obj,'Add Area',select3DRadiusArea(obj));

             case 'Brush Deselection'
                 if ~strcmp(obj.ActiveKey,'space'), obj.Action = 'none'; return; end
                 areaindex = select3DRadiusArea(obj);
                 updateAreaSelectionRadius(obj,'Add Area', areaindex); 
                 if~isempty(areaindex), updateAreaSelection(obj,'Remove Area',areaindex); end
             case 'Brush Dummy'
                 if ~strcmp(obj.ActiveKey,'space'), obj.Action = 'none'; return; end
                 updateAreaSelectionRadius(obj,'Add Area', select3DRadiusArea(obj));
             otherwise
                 return;
         end
         if obj.Record && ~strcmp(obj.Status,'Recording')
            captureFrame(obj);
         end
end

function ScrollWeel(varargin)
         obj = get(varargin{1},'UserData');
         if strcmp(obj.SelectionMode,'brush')&& strcmp(obj.ActiveKey,'space')
            if varargin{2}.VerticalScrollCount < 0
                 obj.SelectionSphere.Radius  = 1.1*obj.SelectionSphere.Radius;
            else
                 obj.SelectionSphere.Radius  = 0.9*obj.SelectionSphere.Radius;
            end
            updateAreaSelectionRadius(obj,'Add Area', select3DRadiusArea(obj));
            if obj.Record
                captureFrame(obj);
            end
            return; 
         end
         if strcmp(obj.SelectionMode,'fill')&& strcmp(obj.ActiveKey,'space')
            if isempty(obj.UserData), return; end
            if varargin{2}.VerticalScrollCount < 0
                 obj.UserData.TH  = 1.1*obj.UserData.TH;
            else
                 obj.UserData.TH   = 0.9*obj.UserData.TH ;
            end
            updateAreaSelection(obj,'Set Area', find(obj.UserData.Distances<=obj.UserData.TH));
            if obj.Record
                captureFrame(obj);
            end
            return; 
         end
         if strcmp(obj.ActiveKey,'space'), return; end
         switch obj.Mode;
             case 'camera'
                 if varargin{2}.VerticalScrollCount < 0
                    q = 1.2;
                    obj.Action = 'zoom in';
                 else
                    q = 0.8;
                    obj.Action = 'zoom out';
                 end
                 zoomCamera(obj,[],q);
             case 'light'
                 if ~obj.SceneLightVisible, return; end
                 xy = [0 0];
                 switch obj.SceneLightRotMode
                     case 'vertical'
                         xy(2) = varargin{2}.VerticalScrollCount*10;
                     case 'horizontal'
                         xy(1) = varargin{2}.VerticalScrollCount*10;
                     otherwise
                         return;
                 end
                 rotateLight(obj,xy);
             otherwise
                 return
         end
         obj.Action = 'none';
         if obj.Record
            captureFrame(obj);
         end
end

function val = PointerOnFigure(fig,scrn_pt)
         set(fig,'Units','pixels')
         pos = get(fig,'Position');
         if (scrn_pt(1)>=pos(1))&&(scrn_pt(2)>=pos(2))&&(scrn_pt(1)-pos(1)<=pos(3))&&(scrn_pt(2)-pos(2)<=pos(4))
             val = true;
         else           
             val = false;
         end

end
