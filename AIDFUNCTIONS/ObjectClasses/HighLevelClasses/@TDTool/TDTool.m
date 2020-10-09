classdef TDTool < superClass
    % this class is a class of 3D viewers
    properties
        % figure object, to display sliders and such
        Handles = [];
        ShapeViewer = [];
        BaseCoeff = [];
        BasePred = [];
        CurrentShape = [];
        CurrentCoeff = [];
        CurrentPred = [];
        Tag = 'TD Morpheus'; 
        UserData = [];
        PCScale = 3;
        MaxPC = [];
        PredScale = 1;
        ActivePC = 0;
        ActivePred = 0;
    end
    properties (Transient = true)
        ShapeModel = [];
        PredModel = [];
    end
    properties (Dependent = true)
        LMModel = [];
    end
    methods %Constructor
        function obj = TDTool(varargin)
          obj = obj@superClass(varargin{:});
          obj.ShapeViewer = viewer3DObj;
          obj.ShapeViewer.Application = obj;
          set(obj.ShapeViewer.Toolbar.landmark_mode,'State','on');
          set(obj.ShapeViewer.Toolbar.colorbar,'State','on');
          obj.ShapeViewer.AxesVisible = true;
          obj.ShapeViewer.AxesGrid = true;
          obj.ShapeViewer.PrintValue = false;
          constructTDTool(obj);
          updateShapeModel(obj);
          set(obj.ShapeViewer.Figure,'CloseRequestFcn', {});
          warning off;
        end % 3D viewer Constructor
    end  
    methods % Special Setting and Getting
        function out = get.Handles(obj)
            out = obj.Handles;
        end
        function out = get.ShapeViewer(obj)
            out = obj.ShapeViewer;
            if ~superClass.isH(out), out = []; end
        end
        function out = get.CurrentShape(obj)
            out = obj.CurrentShape;
            if ~superClass.isH(out), out = []; end
        end
        function out = get.ShapeModel(obj)
            out = obj.ShapeModel;
            if ~superClass.isH(out), out = []; end
        end
        function obj = set.ShapeModel(obj,in)
            obj.ShapeModel = in;
            updateShapeModel(obj);
        end
        function out = get.PredModel(obj)
            out = obj.PredModel;
            if ~superClass.isH(out), out = []; end
        end
        function obj = set.PredModel(obj,in)
            obj.PredModel = in;
            obj.ShapeModel = in.Model;
            updatePredModel(obj);
        end
        function obj = set.PCScale(obj,in)
            obj.PCScale = in;
            updateShapeSlider(obj);
        end
        function obj = set.ActivePC(obj,in)
            obj.ActivePC = in;
            updateShapeSlider(obj);
        end
        function obj = set.PredScale(obj,in)
            obj.PredScale = in;
            updatePredSlider(obj);
        end
        function obj = set.ActivePred(obj,in)
            obj.ActivePred = in;
            updatePredSlider(obj);
        end
        function out = get.LMModel(obj)
                 out = obj.ShapeViewer.LandmarkSelection; 
        end
    end
    methods % Interface functions
        function updateShapeTable(obj,row,value)
            if nargin==1% complete table update
               if isempty(obj.ShapeModel)
                   Data = {};
               else
                   Data = cell(obj.ShapeModel.nrEV,5);
                   for i=1:1:obj.ShapeModel.nrEV
                       if i==obj.ActivePC
                          Data{i,1} = true;
                       else
                          Data{i,1} = false;
                       end
                       Data{i,2} = num2str(i);
                       Data{i,3} = obj.CurrentCoeff(i);
                       Data{i,4} = 0;
                       Data{i,5} = obj.ShapeModel.EigStd(i); 
                   end
               end
            else% single row update
                Data = get(obj.Handles.PC.Table,'Data');
                if ~row==0
                         if row==obj.ActivePred
                                Data{row,1} = true;
                         else
                                Data{row,1} = false;
                         end
                         Data{row,3} = value;
                end
            end   
            set(obj.Handles.PC.Table,'Data',Data);
        end
        function updatePredTable(obj,row,value)
            if nargin==1% complete table update
               if isempty(obj.PredModel)
                   Data = {};
               else
                   Data = cell(obj.PredModel.nrX,5);
                   for i=1:1:obj.PredModel.nrX
                       if i==obj.ActivePred
                          Data{i,1} = true;
                       else
                          Data{i,1} = false;
                       end
                       Data{i,2} = obj.PredModel.XNames{i};
                       Data{i,3} = obj.CurrentPred(i);
                       Data{i,4} = obj.PredModel.XAvg(i);
                       Data{i,5} = obj.PredModel.XStd(i); 
                   end
               end
            else% single row update
                Data = get(obj.Handles.Pred.Table,'Data');
                if ~row==0
                         if row==obj.ActivePred
                                Data{row,1} = true;
                         else
                                Data{row,1} = false;
                         end
                         Data{row,3} = value;
                end
            end   
            set(obj.Handles.Pred.Table,'Data',Data);
        end
        function updateTDTable(obj)
            if isempty(obj.ShapeModel), return; end
            if isempty(obj.ShapeViewer.LandmarkSelection)||isempty(obj.ShapeViewer.LandmarkSelection.Vertices)
                Data = {};       
            else
                nrLM = obj.ShapeViewer.LandmarkSelection.nrV;
                Data = cell(nrLM,4); 
                for i=1:1:nrLM
                    Data{i,1} = obj.ShapeViewer.LandmarkSelection.Vertices(1,i);
                    Data{i,2} = obj.ShapeViewer.LandmarkSelection.Vertices(2,i);
                    Data{i,3} = obj.ShapeViewer.LandmarkSelection.Vertices(3,i);
                    Data{i,4} = obj.ShapeViewer.LandmarkSelection.Value(i);
                end           
            end
            set(obj.Handles.TD.Table,'Data',Data);
        end
        function updateShapeSlider(obj)
            if isempty(obj.ShapeModel), return; end
            if obj.ActivePC == 0;
               set(obj.Handles.PC.Slider,'Enable','off');
            else
               set(obj.Handles.PC.Slider,'Enable','on');
               value = obj.CurrentCoeff(obj.ActivePC);
               mean = 0;
               std = obj.ShapeModel.EigStd(obj.ActivePC);
               set(obj.Handles.PC.Slider,'Min',mean-obj.PCScale*std,'Max',mean+obj.PCScale*std);       
               set(obj.Handles.PC.Slider,'Value',value);
            end
        end
        function updatePredSlider(obj)
            if isempty(obj.PredModel), return; end
            if obj.ActivePred == 0;
               set(obj.Handles.Pred.Slider,'Enable','off');
            else
               set(obj.Handles.Pred.Slider,'Enable','on');
               value = obj.CurrentPred(obj.ActivePred);
               mean = obj.PredModel.XAvg(obj.ActivePred);
               std = obj.PredModel.XStd(obj.ActivePred);
               set(obj.Handles.Pred.Slider,'Min',mean-obj.PredScale*std,'Max',mean+obj.PredScale*std);       
               set(obj.Handles.Pred.Slider,'Value',value);
            end
        end
        function updateShape(obj)
                 scan = getScan(obj.ShapeModel,obj.CurrentCoeff);
                 obj.CurrentShape.Vertices = scan.Vertices;
                 if isempty(scan.TextureMap)
                     obj.CurrentShape.TextureColor = scan.TextureColor;
                 else
                     obj.CurrentShape.TextureMap = clone(scan.TextureMap);
                 end
                 if ~isempty(scan.Value)
                     obj.CurrentShape.Value = scan.Value; %#ok<MCHV3>
                 end
                 %obj.CurrentFace.UserData = scan.UserData;
                 delete(scan);
        end
        function updatePredModel(obj)
            if isempty(obj.PredModel)
                %ActivateAll(obj,'off');
                return;
            end
            obj.BaseCoeff = obj.PredModel.YAvg;
            obj.BasePred = obj.PredModel.XAvg;
            delete(obj.CurrentShape);
            obj.CurrentShape = clone(obj.PredModel.Average);
            obj.CurrentCoeff = obj.PredModel.YAvg;
            obj.CurrentPred = obj.PredModel.XAvg;
            obj.CurrentShape.Axes = obj.ShapeViewer.RenderAxes;
            obj.CurrentShape.Selected = true;
            obj.CurrentShape.Visible = true;
            updateShapeTable(obj);
            updatePredTable(obj);
        end
          function updateShapeModel(obj)
            if isempty(obj.ShapeModel)
                %ActivateAll(obj,'off');
                return;
            end
            obj.BaseCoeff = zeros(1,obj.ShapeModel.nrEV);
            obj.CurrentShape = clone(obj.ShapeModel.Average);
            obj.CurrentCoeff = zeros(1,obj.ShapeModel.nrEV);
            obj.CurrentShape.Axes = obj.ShapeViewer.RenderAxes;
            obj.CurrentShape.Selected = true;
            obj.CurrentShape.Visible = true;
            if strcmp(obj.ShapeModel.Type,'ShapeValuePCA')
               set(obj.ShapeViewer.RenderAxes,'clim',obj.ShapeModel.Range);
            end
            updateShapeTable(obj);
        end
        

    end
    methods % Delete
        function delete(obj)
           if ishandle(obj.Handles.Panel), delete(obj.Handles.Panel);end
           delete@superClass(obj);
           warning on;
        end
    end
end % classdef


%         function ActivatePCPanel(obj,in)
%             set(obj.Handles.PC.Slider,'Enable',in);
%             set(obj.Handles.PC.Popupmenu,'Enable',in);
%             set(obj.Handles.PC.Edit,'Enable',in);
%         end
%         function ActivatePROPPanel(obj,in)
%             set(obj.Handles.PROP.Slider,'Enable',in);
%             set(obj.Handles.PROP.Popupmenu,'Enable',in);
%             set(obj.Handles.PROP.Edit,'Enable',in);
%             set(obj.Handles.PROP.Listbox,'Enable',in);
%             set(obj.Handles.PROP.Popupmenu2,'Enable',in);
%             set(obj.Handles.PROP.PopupmenuScale,'Enable',in);
%             set(obj.Handles.PROP.Button,'Enable',in);
%         end
%         function ActivateCARPanel(obj,in)
%             set(obj.Handles.CAR.Slider,'Enable',in);
%             set(obj.Handles.CAR.Popupmenu,'Enable',in);
%             set(obj.Handles.CAR.Edit,'Enable',in);
%         end
%         function ActivateButton(obj,in)
%             set(obj.Handles.ResetButton,'Enable',in);
%             set(obj.Handles.ExportButton,'Enable',in);
%             set(obj.Handles.ImportButton,'Enable',in);
%             set(obj.Handles.RandomButton,'Enable',in);
%             set(obj.Handles.AverageButton,'Enable',in);
%             set(obj.Handles.BaseButton,'Enable',in);
%         end
%         function ActivateAll(obj,in)
%             ActivatePCPanel(obj,in);
%             ActivatePROPPanel(obj,in)
%             ActivateButton(obj,in);
%             ActivateCARPanel(obj,in)
%         end             