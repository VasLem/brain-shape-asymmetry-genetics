function imageMorphAnalysis(analysis,RefScanV,vref,name,mi,pl,saveresults)
%%      imaging the scan
         RefScan = clone(RefScanV);delete(RefScan.PoseLM);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Vertices = analysis.morphs(2).scan.Vertices;
         if ~isempty(analysis.morphs(2).scan.PoseLM)
            RefScan.PoseLM = clone(analysis.morphs(2).scan.PoseLM);
         end
         RefScan.ColorMode = 'Single';
         set(v.Figure,'Name',[name ' ' pl]);set(v.Figure,'Position',get(vref.Figure,'Position'));
         v.SceneLightPosition = vref.SceneLightPosition;
         drawnow;
         if saveresults
          save(get(v.Figure,'Name'),'RefScan');
          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
          
         end
%%      imaging the norm
         RefScan = clone(RefScanV);delete(RefScan.PoseLM);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Vertices = analysis.morphs(1).scan.Vertices;
         if ~isempty(analysis.morphs(1).scan.PoseLM)
            RefScan.PoseLM = clone(analysis.morphs(1).scan.PoseLM);
         end
         RefScan.ColorMode = 'Single';
         set(v.Figure,'Name',[name ' ' mi]);set(v.Figure,'Position',get(vref.Figure,'Position'));
         v.SceneLightPosition = vref.SceneLightPosition;
         drawnow;
         if saveresults
          save(get(v.Figure,'Name'),'RefScan');
          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
          
         end
%%       imaging the effect map
         RefScan = clone(RefScanV);delete(RefScan.PoseLM);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = analysis.STATS.LocalE(4,:);
         RefScan.ColorMode = 'Indexed';
         set(v.RenderAxes,'clim',[0 max(RefScan.Value)]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name ' Effect']);set(v.Figure,'Position',get(vref.Figure,'Position'));
         v.SceneLightPosition = vref.SceneLightPosition;
         drawnow;
         if saveresults
          save(get(v.Figure,'Name'),'RefScan');
          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
          
         end
%%       imaging the effect-size map
         RefScan = clone(RefScanV);delete(RefScan.PoseLM);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = analysis.STATS.LocalS;
         RefScan.ColorMode = 'Indexed';
         set(v.RenderAxes,'clim',[0 max(RefScan.Value)]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name ' Effect Size']);set(v.Figure,'Position',get(vref.Figure,'Position'));
         v.SceneLightPosition = vref.SceneLightPosition;
         drawnow;
         if saveresults
          save(get(v.Figure,'Name'),'RefScan');
          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
          
         end
%%       imaging the distance map
%          RefScan = clone(RefScanV);
%          v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
%          RefScan.Value = analysis.comparison.Distances;
%          RefScan.ColorMode = 'Indexed';
%          set(v.RenderAxes,'clim',[0 max(RefScan.Value)]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
%          set(v.Figure,'Name',[name ' DistanceMap']);set(v.Figure,'Position',get(vref.Figure,'Position'));
% %          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
%%       X disp
         RefScan = clone(RefScanV);delete(RefScan.PoseLM);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = analysis.comparison.Differences(1,:);
         RefScan.ColorMode = 'Indexed';
         set(v.RenderAxes,'clim',[-1*max(abs(RefScan.Value)) max(abs(RefScan.Value))]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name ' X (Blue Left disp; Red Right disp)']);set(v.Figure,'Position',get(vref.Figure,'Position'));
         drawnow;
         if saveresults
          save(get(v.Figure,'Name'),'RefScan');
          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
          
         end
%%       Y disp
         RefScan = clone(RefScanV);delete(RefScan.PoseLM);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = analysis.comparison.Differences(2,:);
         RefScan.ColorMode = 'Indexed';
         set(v.RenderAxes,'clim',[-1*max(abs(RefScan.Value)) max(abs(RefScan.Value))]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name 'Y (Blue Downward disp; Red Upward disp)']);set(v.Figure,'Position',get(vref.Figure,'Position'));
         drawnow;
         if saveresults
          save(get(v.Figure,'Name'),'RefScan');
          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
          
         end
%%       Z disp
         RefScan = clone(RefScanV);delete(RefScan.PoseLM);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = analysis.comparison.Differences(3,:);
         RefScan.ColorMode = 'Indexed';
         set(v.RenderAxes,'clim',[-1*max(abs(RefScan.Value)) max(abs(RefScan.Value))]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name 'Z (Blue Backwards disp; Red Forward disp)']);set(v.Figure,'Position',get(vref.Figure,'Position'));drawnow;
         if saveresults
          save(get(v.Figure,'Name'),'RefScan');
          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
          
         end    
%%       Normal Disp.
         RefScan = clone(RefScanV);delete(RefScan.PoseLM);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = analysis.comparison.NormalDistances;
         RefScan.ColorMode = 'Indexed';
         set(v.RenderAxes,'clim',[-1*max(abs(RefScan.Value)) max(abs(RefScan.Value))]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name ' Normal (Blue Inward disp; Red Outward disp)']);set(v.Figure,'Position',get(vref.Figure,'Position'));
         drawnow;
         if saveresults
          save(get(v.Figure,'Name'),'RefScan');
          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
          
         end
%%       Area Changes
         RefScan = clone(RefScanV);delete(RefScan.PoseLM);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = analysis.comparison.AreaRatios;
         RefScan.ColorMode = 'Indexed';
         set(v.RenderAxes,'clim',[-1*max(abs(RefScan.Value)) max(abs(RefScan.Value))]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name ' Areal (Blue Decrease; Red Increase )']);set(v.Figure,'Position',get(vref.Figure,'Position'));drawnow;
         if saveresults
          save(get(v.Figure,'Name'),'RefScan');
          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
          
         end      
%%       Curvature Changes
         maxcurv = min(max(analysis.comparison.Curvature1),max(analysis.comparison.Curvature2));
         RefScan = clone(RefScanV);delete(RefScan.PoseLM);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = analysis.comparison.Curvature2;
         RefScan.Vertices = analysis.morphs(2).scan.Vertices;
         RefScan.ColorMode = 'Indexed';
         set(v.RenderAxes,'clim',[0 maxcurv]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name ' Mean Curv. ' pl]);set(v.Figure,'Position',get(vref.Figure,'Position'));
         drawnow;
         if saveresults
          save(get(v.Figure,'Name'),'RefScan');
          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
          
         end
         RefScan = clone(RefScanV);delete(RefScan.PoseLM);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = analysis.comparison.Curvature1;
         RefScan.ColorMode = 'Indexed';
         RefScan.Vertices = analysis.morphs(1).scan.Vertices;
         set(v.RenderAxes,'clim',[0 maxcurv]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name ' Mean Curv. ' mi]);set(v.Figure,'Position',get(vref.Figure,'Position'));
         drawnow;
         if saveresults
          save(get(v.Figure,'Name'),'RefScan');
          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
          
         end
         RefScan = clone(RefScanV);delete(RefScan.PoseLM);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = analysis.comparison.CurvatureRatios;
         RefScan.ColorMode = 'Indexed';
         set(v.RenderAxes,'clim',[-1*max(abs(RefScan.Value)) max(abs(RefScan.Value))]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name ' Mean Curv. (Blue Decrease; Red Increase )']);set(v.Figure,'Position',get(vref.Figure,'Position'));
         drawnow;
         if saveresults
          save(get(v.Figure,'Name'),'RefScan');
          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
          
         end
end