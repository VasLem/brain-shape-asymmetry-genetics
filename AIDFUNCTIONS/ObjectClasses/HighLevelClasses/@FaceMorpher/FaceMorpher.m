classdef FaceMorpher < superClass
    % this class is a class of 3D viewers
    properties
        % figure object, to display sliders and such
        Handles = [];
        Viewer = [];
        BaseCoeff = [];
        CurrentFace = [];
        CurrentCoeff = [];
        CarCoeff = [];
        Tag = 'Face Morpher'; 
        UserData = [];
        PropScale = 2; 
    end
    properties (Transient = true)
        Model = [];
    end
    properties (Dependent = true)
        
    end
    methods %Constructor
        function obj = FaceMorpher(varargin)
          obj = obj@superClass(varargin{:});
          obj.Viewer = viewer3DObj;
          constructFaceMorpher(obj);
          updateModel(obj);
          set(obj.Viewer.Figure,'CloseRequestFcn', {});
        end % 3D viewer Constructor
    end  
    methods % Special Setting and Getting
        function out = get.Handles(obj)
            out = obj.Handles;
            %if ~superClass.isH(out), out = []; end
        end
        function out = get.Viewer(obj)
            out = obj.Viewer;
            if ~superClass.isH(out), out = []; end
        end
        function out = get.CurrentFace(obj)
            out = obj.CurrentFace;
            if ~superClass.isH(out), out = []; end
        end
        function out = get.Model(obj)
            out = obj.Model;
            if ~superClass.isH(out), out = []; end
        end
        function obj = set.Model(obj,in)
            obj.Model = in;
            updateModel(obj);
        end
    end
    methods % Interface functions
        function ActivatePCPanel(obj,in)
            set(obj.Handles.PC.Slider,'Enable',in);
            set(obj.Handles.PC.Popupmenu,'Enable',in);
            set(obj.Handles.PC.Edit,'Enable',in);
        end
        function ActivatePROPPanel(obj,in)
            set(obj.Handles.PROP.Slider,'Enable',in);
            set(obj.Handles.PROP.Popupmenu,'Enable',in);
            set(obj.Handles.PROP.Edit,'Enable',in);
            set(obj.Handles.PROP.Listbox,'Enable',in);
            set(obj.Handles.PROP.Popupmenu2,'Enable',in);
            set(obj.Handles.PROP.PopupmenuScale,'Enable',in);
            set(obj.Handles.PROP.Button,'Enable',in);
        end
        function ActivateCARPanel(obj,in)
            set(obj.Handles.CAR.Slider,'Enable',in);
            set(obj.Handles.CAR.Popupmenu,'Enable',in);
            set(obj.Handles.CAR.Edit,'Enable',in);
        end
        function ActivateButton(obj,in)
            set(obj.Handles.ResetButton,'Enable',in);
            set(obj.Handles.ExportButton,'Enable',in);
            set(obj.Handles.ImportButton,'Enable',in);
            set(obj.Handles.RandomButton,'Enable',in);
            set(obj.Handles.AverageButton,'Enable',in);
            set(obj.Handles.BaseButton,'Enable',in);
        end
        function ActivateAll(obj,in)
            ActivatePCPanel(obj,in);
            ActivatePROPPanel(obj,in)
            ActivateButton(obj,in);
            ActivateCARPanel(obj,in)
        end             
    end
    methods % Delete
        function delete(obj)
           if ishandle(obj.Handles.Panel), delete(obj.Handles.Panel);end
           delete@superClass(obj);
        end
    end
end % classdef
