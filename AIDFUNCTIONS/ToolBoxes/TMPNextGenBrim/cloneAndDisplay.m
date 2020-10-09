function [v,RefScan] = cloneAndDisplay(vref,RefScanV,what,value,color,name,range)
         if nargin < 7, range = []; end
         RefScan = clone(RefScanV);delete(RefScan.PoseLM);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         eval(['RefScan.' what ' = value;']);
         switch color
             case 'Single'
                 RefScan.ColorMode = 'Single';
             case 'Indexed'
                  RefScan.ColorMode = 'Indexed';
                  if isempty(range)
                    set(v.RenderAxes,'clim',[-1*max(abs(RefScan.Value)) max(abs(RefScan.Value))]);
                  else
                    set(v.RenderAxes,'clim',[range(1) range(2)]);
                  end
                  colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
             case 'Texture'
                 RefScan.ColorMode = 'Texture';
         end
         set(v.Figure,'Name',name);set(v.Figure,'Position',get(vref.Figure,'Position'));
         v.SceneLightPosition = vref.SceneLightPosition;
         drawnow;
end