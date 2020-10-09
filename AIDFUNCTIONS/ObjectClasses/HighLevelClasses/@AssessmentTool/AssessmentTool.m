classdef AssessmentTool < superClass
    % this class is a class of 3D viewers
    properties
        % figure object, to display sliders and such
        Handles = [];
        DysmorphViewer = [];
        DistanceViewer = [];
        ThresholdViewer = [];
        VectorViewer = [];
        ScanViewer = [];
        NormViewer = [];
        ScanCurvatureViewer = [];
        NormCurvatureViewer = [];
        CurvatureDiffViewer = [];
        AreaViewer = [];
        NormDisplacementViewer = [];
        AssessmentList = {};
        Active = 0;
        Tag = 'Assessment Tool V1.2: Copyright(c) 2013 P.Claes'; 
        UserData = [];
        Status = 'Ready';
    end
    properties (Transient = true)
        
    end
    properties (Dependent = true)
        ActiveAssessment;
        NrAssessments;
    end
    methods %Constructor
        function obj = AssessmentTool(varargin)
          obj = obj@superClass(varargin{:});
          constructAssessmentTool(obj);
        end
    end  
    methods % Special Setting and Getting
        function out = get.DysmorphViewer(obj)
            out = obj.DysmorphViewer;
            if ~superClass.isH(out), out = []; end
        end
        function obj = set.DysmorphViewer(obj,in)
            obj.DysmorphViewer = in;
            in.Tag = 'Outlier Map';
            set(in.Figure,'CloseRequestFcn', @myCloseFunction);
            in.Application = obj;
            set(in.RenderAxes,'clim',[0 1]);
            colorbar('peer',in.RenderAxes);
            colormap(in.RenderAxes,'hot');
        end
        function out = get.DistanceViewer(obj)
            out = obj.DistanceViewer;
            if ~superClass.isH(out), out = []; end
        end
        function obj = set.DistanceViewer(obj,in)
            obj.DistanceViewer = in;
            in.Tag = 'Distance Map';
            set(in.Figure,'CloseRequestFcn', @myCloseFunction);
            in.Application = obj;
            set(in.RenderAxes,'clim',[0 10]);
            colorbar('peer',in.RenderAxes);
        end
        function out = get.ThresholdViewer(obj)
            out = obj.ThresholdViewer;
            if ~superClass.isH(out), out = []; end
        end
        function obj = set.ThresholdViewer(obj,in)
            obj.ThresholdViewer = in;
            in.Tag = 'Threshold Map';
            set(in.Figure,'CloseRequestFcn', @myCloseFunction);
            in.Application = obj;
            set(in.RenderAxes,'clim',[0 1]);
            colorbar('peer',in.RenderAxes);
            colormap(in.RenderAxes,'winter');
        end
        function out = get.VectorViewer(obj)
            out = obj.VectorViewer;
            if ~superClass.isH(out), out = []; end
        end
        function obj = set.VectorViewer(obj,in)
            obj.VectorViewer = in;
            in.Tag = 'Vector Field';
            set(in.Figure,'CloseRequestFcn', @myCloseFunction);
            in.Application = obj;
            set(in.RenderAxes,'clim',[0 10]);
            colorbar('peer',in.RenderAxes);
        end
        function out = get.ScanViewer(obj)
            out = obj.ScanViewer;
            if ~superClass.isH(out), out = []; end
        end
        function obj = set.ScanViewer(obj,in)
            obj.ScanViewer = in;
            in.Tag = 'Scan';
            set(in.Figure,'CloseRequestFcn', @myCloseFunction);
            in.Application = obj;
        end
        function out = get.NormViewer(obj)
            out = obj.NormViewer;
            if ~superClass.isH(out), out = []; end
        end
        function obj = set.NormViewer(obj,in)
            obj.NormViewer = in;
            in.Tag = 'Norm';
            set(in.Figure,'CloseRequestFcn', @myCloseFunction);
            in.Application = obj;
        end
        function out = get.NormCurvatureViewer(obj)
            out = obj.NormCurvatureViewer;
            if ~superClass.isH(out), out = []; end
        end
        function obj = set.NormCurvatureViewer(obj,in)
            obj.NormCurvatureViewer = in;
            in.Tag = 'Norm Curv.';
            set(in.Figure,'CloseRequestFcn', @myCloseFunction);
            in.Application = obj;
        end
        function out = get.ScanCurvatureViewer(obj)
            out = obj.ScanCurvatureViewer;
            if ~superClass.isH(out), out = []; end
        end
        function obj = set.ScanCurvatureViewer(obj,in)
            obj.ScanCurvatureViewer = in;
            in.Tag = 'Scan Curv.';
            set(in.Figure,'CloseRequestFcn', @myCloseFunction);
            in.Application = obj;
        end
        function out = get.CurvatureDiffViewer(obj)
            out = obj.CurvatureDiffViewer;
            if ~superClass.isH(out), out = []; end
        end
        function obj = set.CurvatureDiffViewer(obj,in)
            obj.CurvatureDiffViewer = in;
            in.Tag = 'Curv. Diff.';
            set(in.Figure,'CloseRequestFcn', @myCloseFunction);
            in.Application = obj;
        end
        function out = get.NormDisplacementViewer(obj)
            out = obj.NormDisplacementViewer;
            if ~superClass.isH(out), out = []; end
        end
        function obj = set.NormDisplacementViewer(obj,in)
            obj.NormDisplacementViewer = in;
            in.Tag = 'Norm. Displ.';
            set(in.Figure,'CloseRequestFcn', @myCloseFunction);
            in.Application = obj;
        end
        function out = get.AreaViewer(obj)
            out = obj.AreaViewer;
            if ~superClass.isH(out), out = []; end
        end
        function obj = set.AreaViewer(obj,in)
            obj.AreaViewer = in;
            in.Tag = 'Area Ratio';
            set(in.Figure,'CloseRequestFcn', @myCloseFunction);
            in.Application = obj;
        end
        function out = get.ActiveAssessment(obj)
            if obj.Active == 0, out = []; return; end
            out = obj.AssessmentList{obj.Active};
            if ~superClass.isH(out), out = []; end
        end
        function obj = set.ActiveAssessment(obj,in)
            if obj.Active == 0, return; end
            obj.AssessmentList{obj.Active} = in;
            updateActiveAssessment(obj);
        end
        function obj = set.Active(obj,in)
            obj.Active = in;
            updateViewers(obj);
            updateActionPanel(obj);
            updateSettingPanel(obj);
                
        end
        function out = get.AssessmentList(obj)
            out = obj.AssessmentList;
        end
        function out = get.NrAssessments(obj)
            out = length(obj.AssessmentList);
        end
        function obj = set.Status(obj,in)
            obj.Status = in;
            switch lower(in)
                case 'ready'
                    %disp('ready');
                    set(obj.Handles.Panel,'Pointer','arrow');
                    drawnow;
                case 'busy'
                    %disp('busy');
                    set(obj.Handles.Panel,'Pointer','watch');
                    drawnow;
            end
        end
    end
    methods % Interface functions
        % update Viewers
        function updateViewers(obj)
            obj.Status = 'busy';
            updateDysmorhViewer(obj);
            updateDistanceViewer(obj);
            updateThresholdViewer(obj);
            updateVectorViewer(obj);
            updateScanViewer(obj);
            updateNormViewer(obj);
            updateScanCurvatureViewer(obj);
            updateNormCurvatureViewer(obj);
            updateCurvatureDiffViewer(obj);
            updateNormDisplacementViewer(obj);
            updateAreaViewer(obj);
            obj.Status = 'ready';           
        end
        function updateDysmorhViewer(obj)
            if isempty(obj.DysmorphViewer), return; end
            if ~obj.DysmorphViewer.Visible, return; end
            deletePatchChildren(obj.DysmorphViewer);
            if obj.Active == 0; return; end
            showDysmorphogram(obj.ActiveAssessment,obj.DysmorphViewer);
        end
        function updateDistanceViewer(obj)
            if isempty(obj.DistanceViewer), return; end
            if ~obj.DistanceViewer.Visible, return; end
            deletePatchChildren(obj.DistanceViewer);
            if obj.Active == 0; return; end
            showDistanceMap(obj.ActiveAssessment,obj.DistanceViewer);
        end
        function updateThresholdViewer(obj)
            if isempty(obj.ThresholdViewer), return; end
            if ~obj.ThresholdViewer.Visible, return; end
            deletePatchChildren(obj.ThresholdViewer);
            if obj.Active == 0; return; end
            showThresholdMap(obj.ActiveAssessment,obj.ThresholdViewer);
        end
        function updateVectorViewer(obj)
            if isempty(obj.VectorViewer), return; end
            if ~obj.VectorViewer.Visible, return; end
            deletePatchChildren(obj.VectorViewer);
            if obj.Active == 0; return; end
            showVectorField(obj.ActiveAssessment,obj.VectorViewer);
        end
        function updateScanViewer(obj)
            if isempty(obj.ScanViewer), return; end
            if ~obj.ScanViewer.Visible, return; end
            hidePatchChildren(obj.ScanViewer);
            if obj.Active == 0; return; end
            showScan(obj.ActiveAssessment,obj.ScanViewer);
        end
        function updateNormViewer(obj)
            if isempty(obj.NormViewer), return; end
            if ~obj.NormViewer.Visible, return; end
            hidePatchChildren(obj.NormViewer);
            if obj.Active == 0; return; end
            showNorm(obj.ActiveAssessment,obj.NormViewer);
        end
        function updateNormCurvatureViewer(obj)
            if isempty(obj.NormCurvatureViewer), return; end
            if ~obj.NormCurvatureViewer.Visible, return; end
            hidePatchChildren(obj.NormCurvatureViewer);
            if obj.Active == 0; return; end
            showNormCurvature(obj.ActiveAssessment,obj.NormCurvatureViewer);
        end
        function updateScanCurvatureViewer(obj)
            if isempty(obj.ScanCurvatureViewer), return; end
            if ~obj.ScanCurvatureViewer.Visible, return; end
            hidePatchChildren(obj.ScanCurvatureViewer);
            if obj.Active == 0; return; end
            showScanCurvature(obj.ActiveAssessment,obj.ScanCurvatureViewer);
        end
        function updateCurvatureDiffViewer(obj)
            if isempty(obj.CurvatureDiffViewer), return; end
            if ~obj.CurvatureDiffViewer.Visible, return; end
            hidePatchChildren(obj.CurvatureDiffViewer);
            if obj.Active == 0; return; end
            showCurvatureDiff(obj.ActiveAssessment,obj.CurvatureDiffViewer);
        end
        function updateNormDisplacementViewer(obj)
            if isempty(obj.NormDisplacementViewer), return; end
            if ~obj.NormDisplacementViewer.Visible, return; end
            hidePatchChildren(obj.NormDisplacementViewer);
            if obj.Active == 0; return; end
            showNormDisplacement(obj.ActiveAssessment,obj.NormDisplacementViewer);
        end
        function updateAreaViewer(obj)
            if isempty(obj.AreaViewer), return; end
            if ~obj.AreaViewer.Visible, return; end
            hidePatchChildren(obj.AreaViewer);
            if obj.Active == 0; return; end
            showArea(obj.ActiveAssessment,obj.AreaViewer);
        end
        % update Panels
        function updateSettingPanel(obj)
            if obj.Active == 0, return;end
            set(obj.Handles.CompensateScaleBox,'Value',obj.ActiveAssessment.CompensateScale);
            if isnan(obj.ActiveAssessment.Significance)
               set(obj.Handles.RobustnessLevelBox,'Value',0);
               set(obj.Handles.RobustnessLevelEdit,'Enable','off');
            else
               set(obj.Handles.RobustnessLevelBox,'Value',1);
               set(obj.Handles.RobustnessLevelEdit,'Enable','on');
            end
            set(obj.Handles.RobustnessLevelEdit,'String',num2str(obj.ActiveAssessment.Significance));
            set(obj.Handles.OutliersOnlyBox,'Value',obj.ActiveAssessment.OutliersOnly);
            switch obj.ActiveAssessment.OutliersBinary
                case true
                    set(obj.Handles.OutliersBinaryBox,'Value',1);
                    set(obj.Handles.OutliersBinaryEdit,'Enable','on');
                    set(obj.Handles.OutliersBinaryEdit,'String',num2str(obj.ActiveAssessment.OutliersThreshold));
                case false
                    set(obj.Handles.OutliersBinaryBox,'Value',0);
                    set(obj.Handles.OutliersBinaryEdit,'Enable','off');
                    set(obj.Handles.OutliersBinaryEdit,'String','nan');
            end
            if numel(obj.ActiveAssessment.Threshold)==1
               set(obj.Handles.LocalBox,'Value',0);
               set(obj.Handles.LocalEdit,'Enable','off');
               set(obj.Handles.GlobalBox,'Value',1);
               set(obj.Handles.GlobalEdit,'Enable','on');
               set(obj.Handles.GlobalEdit,'String',num2str(obj.ActiveAssessment.Threshold));
            else
               set(obj.Handles.LocalBox,'Value',1);
               set(obj.Handles.LocalEdit,'Enable','on');
               set(obj.Handles.GlobalBox,'Value',0);
               set(obj.Handles.GlobalEdit,'Enable','off');
               set(obj.Handles.LocalEdit,'ForegroundColor',[0 0.502 0]);
               setappdata(obj.Handles.LocalEdit,'LocalThresholds',obj.ActiveAssessment.Threshold);
            end
        end
        function updateActionPanel(obj)
            if obj.Active == 0
                set(obj.Handles.ImportScanButton,'Enable','off');
                set(obj.Handles.ImportNormButton,'Enable','off');
                set(obj.Handles.UpdateButton,'Enable','off');                
            else
                set(obj.Handles.ImportScanButton,'Enable','on');
                set(obj.Handles.ImportNormButton,'Enable','on');
                set(obj.Handles.UpdateButton,'Enable','on');               
            end
            
        end
        % Update Misc
        function updateTable(obj,row)
                 obj.Status = 'busy';
                 if nargin ==1
                     if isempty(obj.AssessmentList)
                        Data = {};%cell(1,8);
                     else
                        Data = cell(obj.NrAssessments,7);
                         for i=1:1:obj.NrAssessments
                             if i==obj.Active
                                Data{i,1} = true;
                             else
                                Data{i,1} = false;
                             end
                             Data{i,2} = obj.AssessmentList{i}.Update;
                             Data{i,3} = obj.AssessmentList{i}.Tag;
                             if ~isempty(obj.AssessmentList{i}.Scan),Data{i,4} = obj.AssessmentList{i}.Scan.Tag;end
                             if ~isempty(obj.AssessmentList{i}.Norm),Data{i,5} = obj.AssessmentList{i}.Norm.Tag;end
                             Data{i,6} = obj.AssessmentList{i}.PercDysmorph;
                             Data{i,7} = obj.AssessmentList{i}.PercThresh;
                             Data{i,8} = obj.AssessmentList{i}.RMSE;
                         end                        
                     end
                 else
                     Data = get(obj.Handles.Table,'Data');
                     if ~row==0
                         if row==obj.Active
                                Data{row,1} = true;
                         else
                                Data{row,1} = false;
                         end
                         Data{row,2} = obj.AssessmentList{row}.Update;
                         Data{row,3} = obj.AssessmentList{row}.Tag;
                         if ~isempty(obj.AssessmentList{row}.Scan),Data{row,4} = obj.AssessmentList{row}.Scan.Tag;end
                         if ~isempty(obj.AssessmentList{row}.Norm),Data{row,5} = obj.AssessmentList{row}.Norm.Tag;end
                         Data{row,6} = obj.AssessmentList{row}.PercDysmorph;
                         Data{row,7} = obj.AssessmentList{row}.PercThresh;
                         Data{row,8} = obj.AssessmentList{row}.RMSE;
                         Data{row,9} = obj.AssessmentList{row}.NoiseLevel;
                     end
                 end
                %disp('Updating again');
                set(obj.Handles.Table,'Data',Data);
                obj.Status = 'ready';
        end
        function updateActiveAssessment(obj)
            updateViewers(obj);
            updateTable(obj,obj.Active);
        end
    end
    methods % Delete
        function delete(obj)
           if ishandle(obj.Handles.Panel), delete(obj.Handles.Panel);end
           delete@superClass(obj);
        end
    end
end % classdef

% Application specific viewer closing function
function myCloseFunction(hObject,eventdata)
         viewer = get(hObject,'UserData');
         obj = viewer.Application;
         switch lower(viewer.Tag)
             case 'outlier map'
                 set(obj.Handles.DysmorphogramBox,'Value',0);  
             case 'threshold map'
                 set(obj.Handles.ThresholdMapBox,'Value',0);  
             case 'distance map'
                 set(obj.Handles.DistanceMapBox,'Value',0);  
             case 'vector field'
                 set(obj.Handles.VectorFieldBox,'Value',0);  
             case 'scan'
                 set(obj.Handles.ScanBox,'Value',0);  
             case 'norm'
                 set(obj.Handles.NormBox,'Value',0);
             case 'scan curv.'
                 set(obj.Handles.ScanCurvatureBox,'Value',0);
             case 'norm curv.'
                 set(obj.Handles.NormCurvatureBox,'Value',0);
             case 'curv. diff.'
                 set(obj.Handles.CurvatureDiffBox,'Value',0);
             case 'norm. displ.'
                 set(obj.Handles.NormDisplacementBox,'Value',0);
             case 'area ratio'
                 set(obj.Handles.AreaBox,'Value',0);
         end         
         viewer.Visible = false;
end
