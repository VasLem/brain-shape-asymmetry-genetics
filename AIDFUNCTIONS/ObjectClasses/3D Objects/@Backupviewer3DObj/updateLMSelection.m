function updateLMSelection(obj,action,varargin)
         switch action
             case 'Add LM'
                 obj.LandmarkSelection.Location = [obj.LandmarkSelection.Location varargin{1}];
                 obj.LandmarkSelection.Vi = [obj.LandmarkSelection.Vi varargin{2}];
                 obj.LandmarkSelection.Fi = [obj.LandmarkSelection.Fi varargin{3}];
             case 'Delete LM'
                 index = findClosestLM(obj.LandmarkSelection,varargin{1},10);
                 if isempty(index), return; end
                 full_index = (1:size(obj.LandmarkSelection.Location,2));
                 rest_index = setdiff(full_index,index);
                 obj.LandmarkSelection.Location = obj.LandmarkSelection.Location(:,rest_index);
                 obj.LandmarkSelection.Vi = obj.LandmarkSelection.Vi(:,rest_index);
                 obj.LandmarkSelection.Fi = obj.LandmarkSelection.Fi(:,1:rest_index);
             case 'Delete Last LM'
                 if isempty(obj.LandmarkSelection.Location), return; end% nothing to do
                 if size(obj.LandmarkSelection.Location,2) == 1% patch is going to be invalid
                     updateLMSelection('Delete All',obj); return;
                 end
                 obj.LandmarkSelection.Location = obj.LandmarkSelection.Location(:,1:end-1);
                 obj.LandmarkSelection.Vi = obj.LandmarkSelection.Vi(:,1:end-1);
                 obj.LandmarkSelection.Fi = obj.LandmarkSelection.Fi(:,1:end-1);
             case 'Clear'
                 obj.LandmarkSelection.Location = [];
                 obj.LandmarkSelection.Vi = [];
                 obj.LandmarkSelection.Fi = [];
             case 'Link Pose LM'
                 %if size(obj.LandmarkSelection.Location,2) == 5
                     obj.LandmarkSelection.Tag = 'Pose Landmarks';
                     obj.CurrentMesh.PoseLM = obj.LandmarkSelection;
                     obj.CurrentMesh.PoseLM.SingleColor = [0 0 1];
                     obj.LandmarkSelection = LMObj('Visible',false);
                     obj.LandmarkSelection.SingleColor = [1 0 0];
                     if ~isempty(obj.Parent)
                         handles = guidata(obj.Parent);
                         handles.viewerChange(handles,'Link Pose LM',obj.CurrentMesh);
                     end
                 %else
                 %    msgbox('Exactly 5 Pose Landmarks are required')
                 %end
             case 'Link Custom LM'
%                  obj.LandmarkSelection.Tag = 'Pose Landmarks';
%                  obj.CurrentMesh.PoseLM = obj.LandmarkSelection;
%                  obj.CurrentMesh.PoseLM.SingleColor = [0 0 1];
%                  obj.LandmarkSelection = LMObj('Visible',false);
%                  obj.LandmarkSelection.SingleColor = [1 0 0];
%                  if isempty(obj.Parent), return; end
%                  handles = guidata(obj.Parent);
%                  handles.viewerChange(handles,'Link Custom LM',obj.CurrentMesh);
             otherwise
                 return;
         end
         if ~isempty(obj.LandmarkSelection)
             if size(obj.LandmarkSelection.Location,2) > 0
                obj.LandmarkSelection.Visible = true;
                render(obj.LandmarkSelection);
             else
                obj.LandmarkSelection.Visible = false;
             end
         end
         drawnow;
end