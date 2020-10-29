classdef GeneTool < superClass
    % this class is a class of 3D viewers
    properties
        % figure object, to display sliders and such
        Handles = [];
        ShapeViewer = [];
        BaseCoeff = [];
        BaseC = [];
        BaseX = [];
        CurrentShape = [];
        CurrentCoeff = [];
        CurrentC = [];
        CurrentX = [];
        Tag = 'Janus Morpheus'; 
        UserData = [];
        PCScale = 3;
        MaxPC = [];
        CScale = 3;
        XScale = 4;
        ActivePC = 0;
        ActiveC = 0;
        Curvature = [];
        Area = [];
        Normal = [];
        RenderMode = 'Shape';
        %ActiveX = 0;
    end
    properties (Dependent = true)
        ActiveX;
    end
    properties (Transient = true)
        ShapeModel = [];
        PredModel = [];
    end
    methods %Constructor
        function obj = GeneTool(varargin)
          obj = obj@superClass(varargin{:});
          obj.ShapeViewer = viewer3DObj;
          obj.ShapeViewer.Application = obj;
          set(obj.ShapeViewer.Figure,'CloseRequestFcn', {});
          set(obj.ShapeViewer.Toolbar.light_toggle,'State','on');
          set(obj.ShapeViewer.Toolbar.link_toggle,'State','on');
          constructGeneTool(obj);
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
            %updateShapeModel(obj);
        end
        function out = get.PredModel(obj)
            out = obj.PredModel;
            if ~superClass.isH(out), out = []; end
        end
        function obj = set.PredModel(obj,in)
            reconnectChildren(in);
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
        function obj = set.CScale(obj,in)
            obj.CScale = in;
            updateCSlider(obj);
        end
        function obj = set.ActiveC(obj,in)
            obj.ActiveC = in;
            updateCSlider(obj);
        end
        function obj = set.XScale(obj,in)
            obj.XScale = in;
            updateXSlider(obj);
        end
        function obj = set.ActiveX(obj,in)
            obj.PredModel.ActiveX = in;
            updateXSlider(obj);
        end
        function out = get.ActiveX(obj)
           if isempty(obj.PredModel), out = 0; return; end
           out = obj.PredModel.ActiveX; 
        end
        function obj = set.RenderMode(obj,in)
            switch lower(in)
                case 'shape'
                    obj.CurrentShape.ColorMode = 'Single';
                    obj.CurrentShape.Value = 0;
                    obj.ShapeViewer.BackgroundColor = [0 0 0];
                case 'appearance'
                    if strcmp(obj.ShapeModel.Type,'shapePCA')
                        obj.CurrentShape.ColorMode = 'Single';
                        obj.CurrentShape.Value = 0;
                        obj.ShapeViewer.BackgroundColor = [0 0 0];
                    else
                        obj.CurrentShape.ColorMode = 'Texture';
                        obj.CurrentShape.Value = 0;
                        obj.ShapeViewer.BackgroundColor = [1 1 1];
                    end
                case 'curvature'
                    if isempty(obj.Curvature), updateCurvature(obj); end
                    obj.CurrentShape.Value = obj.Curvature;
                    obj.CurrentShape.ColorMode = 'Indexed';
                    obj.ShapeViewer.BackgroundColor = [1 1 1];
                case 'area'
                    if isempty(obj.Area), updateArea(obj); end
                    obj.CurrentShape.Value = obj.Area;
                    obj.CurrentShape.ColorMode = 'Indexed';
                    obj.ShapeViewer.BackgroundColor = [1 1 1];
                case 'normal'
                    if isempty(obj.Normal), updateNormal(obj); end
                    obj.CurrentShape.Value = obj.Normal;
                    obj.CurrentShape.ColorMode = 'Indexed';
                    obj.ShapeViewer.BackgroundColor = [1 1 1];
                otherwise
                    error('wrong render mode');
            end
            obj.RenderMode = in;
            m = max(abs(obj.CurrentShape.Value));
            if isnan(m), return; end
            if m==0, return; end
            set(obj.ShapeViewer.RenderAxes,'clim',[-1*m m]);
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
                       Data{i,4} = obj.PredModel.YAvg(i);
                       Data{i,5} = obj.ShapeModel.EigStd(i); 
                   end
               end
            else% single row update
                Data = get(obj.Handles.PC.Table,'Data');
                if ~row==0
                         if row==obj.ActivePC
                                Data{row,1} = true;
                         else
                                Data{row,1} = false;
                         end
                         Data{row,3} = value;
                end
            end   
            set(obj.Handles.PC.Table,'Data',Data);
        end
        function updateCTable(obj,row,value)
            if nargin==1% complete table update
               if isempty(obj.PredModel)
                   Data = {};
               else
                   Data = cell(obj.PredModel.nrC,5);
                   for i=1:1:obj.PredModel.nrC
                       if i==obj.ActiveC
                          Data{i,1} = true;
                       else
                          Data{i,1} = false;
                       end
                       Data{i,2} = obj.PredModel.CNames{i};
                       Data{i,3} = obj.CurrentC(i);
                       Data{i,4} = obj.PredModel.CAvg(i);
                       Data{i,5} = obj.PredModel.CStd(i); 
                   end
               end
            else% single row update
                Data = get(obj.Handles.C.Table,'Data');
                if ~row==0
                         if row==obj.ActiveC
                                Data{row,1} = true;
                         else
                                Data{row,1} = false;
                         end
                         Data{row,3} = value;
                end
            end   
            set(obj.Handles.C.Table,'Data',Data);
        end
        function updateXTable(obj,row,value)
            if nargin==1% complete table update
               if isempty(obj.PredModel)
                   Data = {};
               else
                   Data = cell(obj.PredModel.nrX,5);
                   for i=1:1:obj.PredModel.nrX
                       if i==obj.ActiveX
                          Data{i,1} = true;
                       else
                          Data{i,1} = false;
                       end
                       Data{i,2} = obj.PredModel.XNames{i};
                       Data{i,3} = obj.CurrentX(i);
                       Data{i,4} = obj.PredModel.XAvg(i);
                       Data{i,5} = obj.PredModel.XStd(i); 
                   end
               end
            else% single row update
                Data = get(obj.Handles.X.Table,'Data');
                if ~row==0
                         if row==obj.ActiveX
                                Data{row,1} = true;
                         else
                                Data{row,1} = false;
                         end
                         Data{row,3} = value;
                end
            end   
            set(obj.Handles.X.Table,'Data',Data);
        end
        function updateShapeSlider(obj)
            if isempty(obj.ShapeModel), return; end
            if obj.ActivePC == 0;
               set(obj.Handles.PC.Slider,'Enable','off');
            else
               set(obj.Handles.PC.Slider,'Enable','on');
               value = obj.CurrentCoeff(obj.ActivePC);
               mean = obj.PredModel.YAvg(obj.ActivePC);
               std = obj.ShapeModel.EigStd(obj.ActivePC);
               set(obj.Handles.PC.Slider,'Min',mean-obj.PCScale*std,'Max',mean+obj.PCScale*std);       
               set(obj.Handles.PC.Slider,'Value',value);
            end
        end
        function updateCSlider(obj)
            if isempty(obj.PredModel), return; end
            if obj.ActiveC == 0;
               set(obj.Handles.C.Slider,'Enable','off');
            else
               set(obj.Handles.C.Slider,'Enable','on');
               value = obj.CurrentC(obj.ActiveC);
               mean = obj.PredModel.CAvg(obj.ActiveC);
               std = obj.PredModel.CStd(obj.ActiveC);
               set(obj.Handles.C.Slider,'Min',mean-obj.CScale*std,'Max',mean+obj.CScale*std);       
               set(obj.Handles.C.Slider,'Value',value);
            end
        end
        function updateXSlider(obj)
            if isempty(obj.PredModel), return; end
            if obj.ActiveX == 0;
               set(obj.Handles.X.Slider,'Enable','off');
            else
               set(obj.Handles.X.Slider,'Enable','on');
               value = obj.CurrentX(obj.ActiveX);
               mean = obj.PredModel.XAvg(obj.ActiveX);
               std = obj.PredModel.XStd(obj.ActiveX);
               set(obj.Handles.X.Slider,'Min',mean-obj.XScale*std,'Max',mean+obj.XScale*std);       
               set(obj.Handles.X.Slider,'Value',value);
            end
        end
        function updateShape(obj)
                 scan = getScan(obj.ShapeModel,obj.CurrentCoeff);
                 obj.CurrentShape.Vertices = scan.Vertices;
                 if isempty(scan.TextureMap)
                     obj.CurrentShape.TextureColor = scan.TextureColor;
                 else
                     obj.CurrentShape.TextureMap = clone(scan.TextureMap);
                     tmp = obj.CurrentShape.TextureColor;
                 end
                 if ~isempty(scan.Value)
                     obj.CurrentShape.Value = scan.Value; %#ok<MCHV3>
                 end
                 delete(scan);
        end
        function updatePCBar(obj)
                 Data = get(obj.Handles.PC.Table,'Data');
                 Data = cell2mat(Data(:,3:end));
                 Z = (Data(:,1)-Data(:,2))./Data(:,3);
                 bar(obj.Handles.PC.Axes,1:length(Z),Z,1);
                 set(obj.Handles.PC.Axes,'xlim',[0 obj.PredModel.nrY],'ylim',[-2.5 2.5]); 
        end
        function updatePredModel(obj)
            if isempty(obj.PredModel)
                %ActivateAll(obj,'off');
                return;
            end
            obj.BaseCoeff = obj.PredModel.YAvg;
            obj.BaseC = CfromY(obj.PredModel,obj.BaseCoeff);
            obj.BaseX = XfromYC(obj.PredModel,obj.BaseCoeff,obj.BaseC);
            delete(obj.CurrentShape);
            obj.CurrentShape = getScan(obj.ShapeModel,obj.BaseCoeff);
            obj.CurrentShape.SingleColor = [0.8 0.8 0.8];
            obj.CurrentShape.Material = 'Dull';
            obj.CurrentCoeff = obj.BaseCoeff;
            obj.CurrentC = obj.BaseC;
            obj.CurrentX = obj.BaseX;
            obj.CurrentShape.Axes = obj.ShapeViewer.RenderAxes;
            obj.CurrentShape.Selected = true;
            obj.CurrentShape.Visible = true;
            updateShapeTable(obj);
            updateCTable(obj);
            updateXTable(obj);
            updatePCBar(obj);
            if strcmp(obj.PredModel.Model.Type,'appearancePCA'), set(obj.Handles.RenderMode,'Value',2); obj.RenderMode = 'appearance';end
        end
        function updateCurvature(obj)
                [Cmean1,~,Dir1,Dir2]=curvature(obj.CurrentShape,true);
                % getting signed curvature
                Dir3 = cross(Dir1',Dir2');
                angles = vectorAngle(obj.CurrentShape.Gradient,Dir3);
                signs = ones(1,length(angles));
                signs(find(angles<=90)) = -1;
                signedCmean1 = signs.*Cmean1';
                % smoothing of curvatures
                obj.Curvature = smoothFunction(obj.CurrentShape,signedCmean1',2,'functiondistance');
                %obj.RenderMode = obj.RenderMode;
                
        end
        function updateArea(obj)             
                baseshape = getScan(obj.ShapeModel,obj.BaseCoeff);
                areas2 = localAreas(obj.CurrentShape);
                areas1 = localAreas(baseshape);
                obj.Area = -1*log(areas1./areas2);
                %obj.RenderMode = obj.RenderMode;
        end
        function updateNormal(obj)
                 baseshape = getScan(obj.ShapeModel,obj.BaseCoeff);
                 out = vNormalDistances(baseshape,obj.CurrentShape);
                 if ~isempty(find(isnan(out)))
                    obj.Normal = zeros(1,baseshape.nrV);
                 else
                    obj.Normal = out;
                 end
                 %obj.RenderMode = obj.RenderMode;
        end
        function updateFC(obj)
           setptr(obj.Handles.Panel,'watch'); drawnow;
           updateCurvature(obj);
           updateArea(obj);
           updateNormal(obj);
           obj.RenderMode = obj.RenderMode;
           setptr(obj.Handles.Panel,'arrow'); drawnow;
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