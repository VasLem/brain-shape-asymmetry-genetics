function constructViewerToolbar(obj)
%% CREATION

load viewertoolbarimagesOBJ.mat;
obj.Toolbar.h = uitoolbar(obj.Figure);
obj.Toolbar.open = uipushtool(obj.Toolbar.h,'Cdata',toolbarimages.open,...
                'TooltipString','Open 3D',...
                'Separator','off',...
                'ClickedCallback',@open_Callback,...
                'UserData',obj);
obj.Toolbar.save = uipushtool(obj.Toolbar.h,'Cdata',toolbarimages.save,...
                'TooltipString','Save 3D',...
                'Separator','off',...
                'ClickedCallback',@save_Callback,...
                'UserData',obj);
obj.Toolbar.cam_mode = uitoggletool(obj.Toolbar.h,'Cdata',toolbarimages.orbit,...
                'TooltipString','Toggle Camera Mode',...
                'Separator','on',...
                'OnCallback',@camera_mode_on_Callback,...
                'OffCallback',@camera_mode_off_Callback,...
                'UserData',obj);
obj.Toolbar.ortho = uipushtool(obj.Toolbar.h,'Cdata',toolbarimages.ortho,...
                'TooltipString','Orthographic Projection',...
                'Separator','off',...
                'ClickedCallback',@camera_ortho_Callback,...
                'UserData',obj);
obj.Toolbar.persp = uipushtool(obj.Toolbar.h,'Cdata',toolbarimages.perspective,...
                'TooltipString','Perspective Projection',...
                'Separator','off',...
                'ClickedCallback',@camera_persp_Callback,...
                'UserData',obj);            
obj.Toolbar.light_mode = uitoggletool(obj.Toolbar.h,'Cdata',toolbarimages.orbitlight,...
                'TooltipString','Toggle Light Mode',...
                'Separator','on',...
                'OnCallback',@light_mode_on_Callback,...
                'OffCallback',@light_mode_off_Callback,...
                'UserData',obj);
obj.Toolbar.light_toggle = uitoggletool(obj.Toolbar.h,'Cdata',toolbarimages.light,...
                'TooltipString','Toggle Light',...
                'Separator','off',...
                'OnCallback',@light_on_Callback,...
                'OffCallback',@light_off_Callback,...
                'UserData',obj);
obj.Toolbar.link_toggle = uitoggletool(obj.Toolbar.h,'Cdata',toolbarimages.link,...
                'TooltipString','Link Light',...
                'Separator','off',...
                'OnCallback',@link_on_Callback,...
                'OffCallback',@link_off_Callback,...
                'UserData',obj);
obj.Toolbar.full_mode = uitoggletool(obj.Toolbar.h,'Cdata',toolbarimages.objectcube,...
                'TooltipString','Full Scan Mode',...
                'Separator','on',...
                'OnCallback',@full_mode_on_Callback,...
                'OffCallback',@full_mode_off_Callback,...
                'UserData',obj);            
obj.Toolbar.landmark_mode = uitoggletool(obj.Toolbar.h,'Cdata',toolbarimages.LMindicator,...
                'TooltipString','Landmark Mode',...
                'Separator','off',...
                'OnCallback',@landmark_mode_on_Callback,...
                'OffCallback',@landmark_mode_off_Callback,...
                'UserData',obj);
obj.Toolbar.area_mode = uitoggletool(obj.Toolbar.h,'Cdata',toolbarimages.brush2,...
                'TooltipString','Area Mode',...
                'Separator','off',...
                'OnCallback',@area_mode_on_Callback,...
                'OffCallback',@area_mode_off_Callback,...
                'UserData',obj);
obj.Toolbar.brush_mode = uitoggletool(obj.Toolbar.h,'Cdata',toolbarimages.brush,...
                'TooltipString','Brush Mode',...
                'Separator','off',...
                'OnCallback',@brush_mode_on_Callback,...
                'OffCallback',@brush_mode_off_Callback,...
                'UserData',obj);
obj.Toolbar.fill_mode = uitoggletool(obj.Toolbar.h,'Cdata',toolbarimages.fill,...
                'TooltipString','Fill Mode',...
                'Separator','off',...
                'OnCallback',@fill_mode_on_Callback,...
                'OffCallback',@fill_mode_off_Callback,...
                'UserData',obj);            
obj.Toolbar.colorbar = uitoggletool(obj.Toolbar.h,'Cdata',toolbarimages.colorbar,...
                'TooltipString','Toggle Colorbar',...
                'Separator','on',...
                'OnCallback',@colorbar_on_Callback,...
                'OffCallback',@colorbar_off_Callback,...
                'UserData',obj);
obj.Toolbar.colorbareditor = uipushtool(obj.Toolbar.h,'Cdata',toolbarimages.colorbareditor,...
                'TooltipString','Edit Colorbar',...
                'Separator','off',...
                'ClickedCallback',@colorbar_editor_Callback,...
                'UserData',[]);
obj.Toolbar.imagecapture = uipushtool(obj.Toolbar.h,'Cdata',toolbarimages.camera,...
                'TooltipString','Edit Colorbar',...
                'Separator','on',...
                'ClickedCallback',@imageCaptureCallback,...
                'UserData',obj);
obj.Toolbar.video = uitoggletool(obj.Toolbar.h,'Cdata',toolbarimages.video,...
                'TooltipString','Toggle Colorbar',...
                'Separator','off',...
                'OnCallback',@recordOnCallback,...
                'OffCallback',@recordOffCallback,...
                'UserData',obj);
           
            

end

%% CALLBACKS
function open_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         scan = meshObj.import;
         if isempty(scan), return; end
         if iscell(scan)
             for i=1:1:length(scan)
                scan{i}.Axes = obj.RenderAxes;
                scan{i}.Visible = true;
                scan{i}.Selected = true; 
             end
         else
             scan.Axes = obj.RenderAxes;
             scan.Visible = true;
             scan.Selected = true;
         end                        
end

function save_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if isempty(obj.MeshChildren), return; end
         for i=1:1:length(obj.MeshChildren)
             if obj.MeshChildren{i}.Selected
                 export(obj.MeshChildren{i});
             end
         end                     
end

function camera_mode_on_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         obj.Mode = 'camera';                   
end

function camera_mode_off_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         currentmode = obj.Mode;
         switch currentmode
             case 'camera'
                 obj.Mode = 'none';
             otherwise
         end         
end

function camera_ortho_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         obj.CamProjection = 'orthographic';
end

function camera_persp_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         obj.CamProjection = 'perspective';
end

function light_mode_on_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         obj.Mode = 'light';              
end

function light_mode_off_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         currentmode = obj.Mode;
         switch currentmode
             case 'light'
                 obj.Mode = 'none';
             otherwise
         end        
end

function light_on_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         %disp('turning light on');
         obj.SceneLightVisible = true;
end

function light_off_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         obj.SceneLightVisible = false;
end

function link_on_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         obj.SceneLightLinked = true;
end

function link_off_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         obj.SceneLightLinked = false;
end

function full_mode_on_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         obj.SelectionMode = 'full';                    
end

function full_mode_off_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         currentmode = obj.SelectionMode;         
         switch currentmode
             case 'landmark'
             case 'area'
             case 'fill'
             case 'brush'
             case 'full'
                 obj.SelectionMode = 'none';
             otherwise
                 set(obj.Figure,'renderer',obj.Renderer);
         end        
         
end

function landmark_mode_on_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         set(obj.Figure,'renderer','zbuffer');
         camproj(obj.RenderAxes,'orthographic');
         obj.SelectionMode = 'landmark';
         obj.LandmarkSelection = LMObj('Visible',false);
         obj.UserData = [];
         obj.LandmarkSelection.SingleColor = [1 0 0];                
         
end

function landmark_mode_off_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         currentmode = obj.SelectionMode;
         if validChild(obj,obj.LandmarkSelection), delete(obj.LandmarkSelection); end
         obj.LandmarkSelection = [];
         switch currentmode
             case 'landmark'
                 obj.SelectionMode = 'none';
                 set(obj.Figure,'renderer',obj.Renderer);
             case 'area'
                 setTitle(obj,''); 
             case 'fill'
                 setTitle(obj,''); 
             case 'brush'
                 setTitle(obj,''); 
             case 'full'
                 setTitle(obj,''); 
             otherwise
                 set(obj.Figure,'renderer',obj.Renderer);
         end                
end

function area_mode_on_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         set(obj.Figure,'renderer','zbuffer');
         camproj(obj.RenderAxes,'orthographic');
         obj.SelectionMode = 'area';
         obj.UserData = [];
         linkAreaSelection(obj);                        
end

function area_mode_off_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         currentmode = obj.SelectionMode;
         switch currentmode
             case 'landmark'
                  if validChild(obj,obj.AreaSelection), delete(obj.AreaSelection); end
                  obj.AreaSelection = [];
             case 'area'
                 if validChild(obj,obj.AreaSelection), delete(obj.AreaSelection); end
                 obj.AreaSelection = [];
                 obj.SelectionMode = 'none';
                 set(obj.Figure,'renderer',obj.Renderer);
             case 'fill'
             case 'brush'
             case 'full'
                 if validChild(obj,obj.AreaSelection), delete(obj.AreaSelection); end
                 obj.AreaSelection = [];
                 set(obj.Figure,'renderer',obj.Renderer);
             otherwise
                 set(obj.Figure,'renderer',obj.Renderer);
                 if validChild(obj,obj.AreaSelection), delete(obj.AreaSelection); end
                 obj.AreaSelection = [];
         end        
         
end

function fill_mode_on_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         set(obj.Figure,'renderer','zbuffer');
         camproj(obj.RenderAxes,'orthographic');
         obj.SelectionMode = 'fill';
         obj.UserData = [];
         linkAreaSelection(obj);                 
         
end

function fill_mode_off_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         currentmode = obj.SelectionMode;        
         switch currentmode
             case 'landmark'
                 if validChild(obj,obj.AreaSelection), delete(obj.AreaSelection); end
                 obj.AreaSelection = [];
             case 'area'
             case 'fill'
                 if validChild(obj,obj.AreaSelection), delete(obj.AreaSelection); end
                 obj.AreaSelection = [];
                 obj.SelectionMode = 'none';
                 set(obj.Figure,'renderer',obj.Renderer);
             case 'brush'
             case 'full'
                 if validChild(obj,obj.AreaSelection), delete(obj.AreaSelection); end
                 obj.AreaSelection = [];
                 set(obj.Figure,'renderer',obj.Renderer);
             otherwise
                 if validChild(obj,obj.AreaSelection), delete(obj.AreaSelection); end
                 obj.AreaSelection = [];
                 set(obj.Figure,'renderer',obj.Renderer);
         end         
         
end

function brush_mode_on_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         set(obj.Figure,'renderer','zbuffer');
         camproj(obj.RenderAxes,'orthographic');
         obj.SelectionMode = 'brush';
         linkAreaSelection(obj);
         obj.UserData = [];
         linkAreaSelectionRadius(obj);           
         
end

function brush_mode_off_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         currentmode = obj.SelectionMode;
         switch currentmode
             case 'landmark'
                  if validChild(obj,obj.AreaSelection), delete(obj.AreaSelection); end
                  obj.AreaSelection = [];
                  if validChild(obj,obj.AreaSelectionRadius), delete(obj.AreaSelectionRadius); end
                  obj.AreaSelectionRadius = [];
             case {'area', 'fill'}
                  if validChild(obj,obj.AreaSelectionRadius), delete(obj.AreaSelectionRadius); end
                  obj.AreaSelectionRadius = [];
             case 'brush'
                 if validChild(obj,obj.AreaSelection), delete(obj.AreaSelection); end
                 obj.AreaSelection = [];
                 if validChild(obj,obj.AreaSelectionRadius), delete(obj.AreaSelectionRadius); end
                 obj.AreaSelectionRadius = [];
                 obj.SelectionMode = 'none';
                 set(obj.Figure,'renderer',obj.Renderer);
             case 'full'
                 if validChild(obj,obj.AreaSelection), delete(obj.AreaSelection); end
                 obj.AreaSelection = [];
                 if validChild(obj,obj.AreaSelectionRadius), delete(obj.AreaSelectionRadius); end
                 obj.AreaSelectionRadius = [];
                 set(obj.Figure,'renderer',obj.Renderer);
             otherwise
                 set(obj.Figure,'renderer',obj.Renderer);
                 if validChild(obj,obj.AreaSelection), delete(obj.AreaSelection); end
                 obj.AreaSelection = [];
                 if validChild(obj,obj.AreaSelectionRadius), delete(obj.AreaSelectionRadius); end
                 obj.AreaSelectionRadius = [];
         end        
         
end

function colorbar_on_Callback(hObject,eventdata)
           obj = get(hObject,'UserData');
           colorbar('peer',obj.RenderAxes,'LOCATION','EastOutside');
end

function colorbar_off_Callback(hObject,eventdata)
           obj = get(hObject,'UserData');
           axes(obj.RenderAxes);
           colorbar('off');
end

function colorbar_editor_Callback(hObject,eventdata)
         colormapeditor;
end

function imageCaptureCallback(hObject,eventdata)
         obj = get(hObject,'UserData');
         im = captureImage(obj);
         pause(0.2);
         [filename, pathname, filterindex] = uiputfile( ...
                                             {'*.tiff','Images (*.tiff)'; ...
                                              '*.png','Images (*.png)'; ...
                                              '*.tiff','Images (*.tiff)'; ...
                                              '*.bmp','Images (*.bmp)'; ...
                                              '*.jpg','Images (*.jpg)'; ...
                                              '*.gif','Images (*.gif)'},...
                                              'Save as', get(obj.Figure,'Name'));
          if filename == 0; return; end
          cd(pathname);             
          switch filterindex
              case 1
                  imwrite(im,filename,'tiff','Compression','none','Resolution',600);
                  %saveas(obj.Figure,filename,'eps');
                  %saveas(obj.Figure,filename,'ai');
                  %print(obj.Figure,'-depsc','-r600',filename(1:end-5));
                  %saveas(obj.Figure,filename(1:end-5),'fig');
                  %disp('saved');
              case 2
                  imwrite(im,filename,'png');
              case 3
                  imwrite(im,filename,'bmp');
              case 4
                  imwrite(im,filename,'jpeg');
              case 5
                  imwrite(im,filename,'gif');
          end
end

function recordOnCallback(hObject,eventdata)
         obj = get(hObject,'UserData');
         [filename, pathname, filterindex] = uiputfile( ...
                                             {'*.avi','Avi files (*.avi)'}, ...
                                              'Save as', 'Movie'); %#ok<NASGU>
         if filename == 0, set(obj.Toolbar.video,'State','off'); return; end
         obj.FrameCounter = 0;
         cd(pathname);
         obj.MovieFile = avifile(filename);
         obj.MovieFile.fps = 6;
         obj.MovieFile.compression = 'none';
         set(obj.Figure,'Resize','off');
         obj.Record = true;
         captureFrame(obj);
end

function recordOffCallback(hObject,eventdata)
         obj = get(hObject,'UserData');
         obj.Record = false;
         set(obj.Figure,'Resize','on');
         if isempty(obj.MovieFile), return; end
         
         obj.MovieFile = close(obj.MovieFile);
         obj.MovieFile = [];
         obj.FrameCounter = 0;
end
