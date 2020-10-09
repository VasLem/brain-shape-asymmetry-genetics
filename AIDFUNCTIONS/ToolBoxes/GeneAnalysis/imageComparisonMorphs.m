function imageComparisonMorphs(comp,RefScanV,vref,name,mi,pl)
%%      imaging the scan
         RefScan = clone(RefScanV);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Vertices = comp.morphs(2).scan.Vertices;
         RefScan.ColorMode = 'Single';
         set(v.Figure,'Name',[name ' ' pl ]);set(v.Figure,'Position',get(vref.Figure,'Position'));
         v.SceneLightPosition = vref.SceneLightPosition;
         drawnow;
%          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
%%      imaging the norm
         RefScan = clone(RefScanV);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Vertices = comp.morphs(1).scan.Vertices;
         RefScan.ColorMode = 'Single';
         set(v.Figure,'Name',[name ' ' mi]);set(v.Figure,'Position',get(vref.Figure,'Position'));
         v.SceneLightPosition = vref.SceneLightPosition;
         drawnow;
%          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
%%       imaging the distance map
         RefScan = clone(RefScanV);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = comp.comparison.Distances;
         RefScan.ColorMode = 'Indexed';
         set(v.RenderAxes,'clim',[0 max(RefScan.Value)]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name ' DistanceMap']);set(v.Figure,'Position',get(vref.Figure,'Position'));
%          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
%%       X disp
         RefScan = clone(RefScanV);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = comp.comparison.Differences(1,:);
         RefScan.ColorMode = 'Indexed';
         set(v.RenderAxes,'clim',[-1*max(abs(RefScan.Value)) max(abs(RefScan.Value))]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name ' X (Blue Left disp; Red Right disp)']);set(v.Figure,'Position',get(vref.Figure,'Position'));
%          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
%%       Y disp
         RefScan = clone(RefScanV);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = comp.comparison.Differences(2,:);
         RefScan.ColorMode = 'Indexed';
         set(v.RenderAxes,'clim',[-1*max(abs(RefScan.Value)) max(abs(RefScan.Value))]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name 'Y (Blue Downward disp; Red Upward disp)']);set(v.Figure,'Position',get(vref.Figure,'Position'));
%          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);         
%%       Z disp
         RefScan = clone(RefScanV);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = comp.comparison.Differences(3,:);
         RefScan.ColorMode = 'Indexed';
         set(v.RenderAxes,'clim',[-1*max(abs(RefScan.Value)) max(abs(RefScan.Value))]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name 'Z (Blue Backwards disp; Red Forward disp)']);set(v.Figure,'Position',get(vref.Figure,'Position'));
%          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);    
%%       Normal Disp.
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = comp.comparison.NormalDistances;
         RefScan.ColorMode = 'Indexed';
         set(v.RenderAxes,'clim',[-1*max(abs(RefScan.Value)) max(abs(RefScan.Value))]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name ' Normal (Blue Inward disp; Red Outward disp)']);set(v.Figure,'Position',get(vref.Figure,'Position'));
%          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);         
%%       Area Changes
         RefScan = clone(RefScanV);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = comp.comparison.AreaRatios;
         RefScan.ColorMode = 'Indexed';
         set(v.RenderAxes,'clim',[-1*max(abs(RefScan.Value)) max(abs(RefScan.Value))]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name ' Areal (Blue Decrease; Red Increase )']);set(v.Figure,'Position',get(vref.Figure,'Position'));
%          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);         
%%       Curvature Changes
         maxcurv = min(max(comp.comparison.Curvature1),max(comp.comparison.Curvature2));
         RefScan = clone(RefScanV);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = comp.comparison.Curvature2;
         RefScan.Vertices = comp.morphs(2).scan.Vertices;
         RefScan.ColorMode = 'Indexed';
         set(v.RenderAxes,'clim',[0 maxcurv]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name ' Mean Curv. ' pl]);set(v.Figure,'Position',get(vref.Figure,'Position'));
%          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);         
         RefScan = clone(RefScanV);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = comp.comparison.Curvature1;
         RefScan.ColorMode = 'Indexed';
         RefScan.Vertices = comp.morphs(1).scan.Vertices;
         set(v.RenderAxes,'clim',[0 maxcurv]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name ' Mean Curv. ' mi]);set(v.Figure,'Position',get(vref.Figure,'Position'));
%          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);         
         RefScan = clone(RefScanV);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = comp.comparison.CurvatureRatios;
         RefScan.ColorMode = 'Indexed';
         set(v.RenderAxes,'clim',[-1*max(abs(RefScan.Value)) max(abs(RefScan.Value))]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name ' Mean Curv. (Blue Decrease; Red Increase )']);set(v.Figure,'Position',get(vref.Figure,'Position'));
%          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);         
end