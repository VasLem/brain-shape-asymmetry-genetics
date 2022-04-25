function updateLMSelection(obj,action,varargin)
         switch action
             case 'Add LM'
                 obj.LandmarkSelection.Location = [obj.LandmarkSelection.Location varargin{1}];
                 obj.LandmarkSelection.Vi = [obj.LandmarkSelection.Vi varargin{2}];
                 obj.LandmarkSelection.Fi = [obj.LandmarkSelection.Fi varargin{3}];
                 cart2bary(obj.LandmarkSelection,obj.CurrentMesh);
                 if isempty(obj.CurrentMesh.IndexedColor),
                    obj.LandmarkSelection.Value = [obj.LandmarkSelection.Value nan];
                 else
                    baryValue(obj.LandmarkSelection,obj.CurrentMesh);
                    if obj.PrintValue
                       setTitle(obj,['Value = ' num2str(obj.LandmarkSelection.Value(end))]);
                    end
                 end
                 if ~isempty(obj.Application), ViewerAction(obj.Application,'LM Updated'); end
             case 'Delete LM'
                 index = findClosestLM(obj.LandmarkSelection,varargin{1},10);
                 if isempty(index), return; end
                 full_index = (1:size(obj.LandmarkSelection.Location,2));
                 rest_index = setdiff(full_index,index);
                 obj.LandmarkSelection.Location = obj.LandmarkSelection.Location(:,rest_index);
                 obj.LandmarkSelection.Vi = obj.LandmarkSelection.Vi(:,rest_index);
                 obj.LandmarkSelection.Fi = obj.LandmarkSelection.Fi(:,rest_index);
                 cart2bary(obj.LandmarkSelection,obj.CurrentMesh);
                 if ~isempty(obj.LandmarkSelection.Value)
                    %obj.LandmarkSelection.Value = obj.LandmarkSelection.Value(:,1:rest_index);
                    baryValue(obj.LandmarkSelection,obj.CurrentMesh);
                 end
                 setTitle(obj,'');
                 if ~isempty(obj.Application), ViewerAction(obj.Application,'LM Updated'); end
             case 'Delete Last LM'
                 if isempty(obj.LandmarkSelection.Location), return; end% nothing to do
                 if size(obj.LandmarkSelection.Location,2) == 1% patch is going to be invalid
                     updateLMSelection('Delete All',obj); return;
                 end
                 obj.LandmarkSelection.Location = obj.LandmarkSelection.Location(:,1:end-1);
                 obj.LandmarkSelection.Vi = obj.LandmarkSelection.Vi(:,1:end-1);
                 obj.LandmarkSelection.Fi = obj.LandmarkSelection.Fi(:,1:end-1);
                 cart2bary(obj.LandmarkSelection,obj.CurrentMesh);
                 if ~isempty(obj.LandmarkSelection.Value)
                    %obj.LandmarkSelection.Value = obj.LandmarkSelection.Value(:,1:end-1);
                    baryValue(obj.LandmarkSelection,obj.CurrentMesh);
                 end
                 setTitle(obj,'');
                 if ~isempty(obj.Application), ViewerAction(obj.Application,'LM Updated'); end
             case 'Clear'
                 obj.LandmarkSelection.Location = [];
                 obj.LandmarkSelection.Vi = [];
                 obj.LandmarkSelection.Fi = [];
                 obj.LandmarkSelection.Value = [];
                 obj.LandmarkSelection.Bary = [];
                 setTitle(obj,'');
                 if ~isempty(obj.Application), ViewerAction(obj.Application,'LM Updated'); end
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
                   obj.LandmarkSelection = clone(obj.CurrentMesh.PoseLM);
                   obj.LandmarkSelection
                   delete(obj.CurrentMesh.PoseLM);
                   if isempty(obj.LandmarkSelection.Value)
                       obj.LandmarkSelection.Value = nan*ones(1,obj.LandmarkSelection.nrV);
                   end
                   %obj.LandmarkSelection = LMObj('Visible',false);
                   obj.LandmarkSelection.SingleColor = [1 0 0];
%                  obj.LandmarkSelection.Tag = 'Pose Landmarks';
%                  obj.CurrentMesh.PoseLM = obj.LandmarkSelection;
%                  obj.CurrentMesh.PoseLM.SingleColor = [0 0 1];
%                  obj.LandmarkSelection = LMObj('Visible',false);
%                  obj.LandmarkSelection.SingleColor = [1 0 0];
%                  if isempty(obj.Parent), return; end
%                  handles = guidata(obj.Parent);
%                  handles.viewerChange(handles,'Link Custom LM',obj.CurrentMesh);
             case 'Show ID'
                 index = findClosestLM(obj.LandmarkSelection,varargin{1},10);
                 if isempty(index), return; end
                 setTitle(obj,['ID = ' num2str(index)]);
             case 'Updated Shape'
                 if ~isempty(obj.LandmarkSelection.Vertices)
                     bary2cart(obj.LandmarkSelection,obj.CurrentMesh);
                     if ~isempty(obj.CurrentMesh.IndexedColor),
                        baryValue(obj.LandmarkSelection,obj.CurrentMesh);
                        if obj.PrintValue
                           setTitle(obj,['Value = ' num2str(obj.LandmarkSelection.Value(end))]);
                        end
                     end
                     if ~isempty(obj.Application), ViewerAction(obj.Application,'LM Updated'); end
                 end
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