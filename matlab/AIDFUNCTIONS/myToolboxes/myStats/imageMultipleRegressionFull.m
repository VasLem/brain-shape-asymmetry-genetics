function imageMultipleRegressionFull(Result,RefScan,vref,name)
%% Multiple effect
try
    RefScan = clone(RefScan);
    v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
    RefScan.Value = double((Result.MT.pLocalR2<=0.001))';colormap(v.RenderAxes,'summer');
    set(v.Figure,'Name',[name ' Mult Sign']);set(v.Figure,'Position',get(vref.Figure,'Position'));
    colormap(v.RenderAxes,'jet');
    im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
catch
end
RefScan = clone(RefScan);
v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
RefScan.Value = Result.MT.LocalR2;set(v.Figure,'Name',[name ' Mult R2']);set(v.Figure,'Position',get(vref.Figure,'Position'));
set(v.RenderAxes,'clim',[0 max(RefScan.Value)]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
%set(v.RenderAxes,'clim',[0 0.45]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
colormap(v.RenderAxes,'jet');
im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
end