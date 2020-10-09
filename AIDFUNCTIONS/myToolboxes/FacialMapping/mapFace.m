function [out] = mapFace(Face,Template,RMapper,NRMapper,display)
     %warning off;
     % Step 0, initialize
      out = clone(Template);% make a copy of the Template
      out.UserData = [];
      Face = clone(Face);
      if nargin<5, display = false;end
     % Init visualization if requested 
      if display % Setting up viewers and then press enter to continue;
              v = viewer(out);
              Face.RenderAxes = v.RenderAxes;
              Face.Visible = true;
              Face.ViewMode = 'Wireframe';
              Face.ColorMode = 'single';
              Face.SingleColor = [0.8 0.8 0.8];
              out.ViewMode = 'Wireframe';
              out.SingleColor = [0.1 0.5 0.9];
              str = input('Are you ready? Y/N [Y]: ','s');
              if isempty(str), str = 'y'; end
              switch lower(str)
                  case 'y'
                      disp('Started...');
                  case 'n'
                      return;
              end
      end  
     % Step 1, landmark based initialisation (crude rigid registration)
      if isfield(Template.UserData,'LandmarkSelection')
          if display, disp('STEP1 Crude rigid registration');end
          [~,~,transform] = procrustes(Face.UserData.LandmarkSelection.Vertices,Template.UserData.LandmarkSelection.Vertices,'Scaling',true,'Reflection',false);
          out.Vertices = transform.b*out.Vertices*transform.T + repmat(transform.c(1,:),out.nVertices,1);
          drawnow;
      end
     % Step 2, ICP based fine initialisation (fine rigid registration)
      if ~isempty(RMapper)
          if display, disp('STEP2 Fine rigid registration');end
          RMapper.TargetShape = Face;
          RMapper.FloatingShape = out;
          RMapper.Display = false;
          map(RMapper);
          out = clone(RMapper.FloatingShape);
      end
     % Step 3, non-rigid ICP based mapping (fine non-rigid registration)
     if ~isempty(NRMapper)
          if display, disp('STEP3 Fine non rigid registration');end
          NRMapper.TargetShape = Face;
          NRMapper.FloatingShape = out;
          NRMapper.Display = false;
          map(NRMapper);
          out = clone(NRMapper.FloatingShape);
     end
end