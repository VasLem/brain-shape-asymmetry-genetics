function imageMultipleAndPartialRegression(Result,RefScan,vref,name)
%% Multiple effect
try
    RefScan = clone(RefScan);
    v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
    RefScan.Value = double((Result.MT.pLocalR2<=0.001))';colormap(v.RenderAxes,'summer');
    set(v.Figure,'Name',[name ' Mult Sign']);set(v.Figure,'Position',get(vref.Figure,'Position'));
    colormap(v.RenderAxes,'gray');
    im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
catch
end
RefScan = clone(RefScan);
v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
RefScan.Value = Result.MT.LocalR2;set(v.Figure,'Name',[name ' Mult R2']);set(v.Figure,'Position',get(vref.Figure,'Position'));
%set(v.RenderAxes,'clim',[0 max(RefScan.Value)]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
set(v.RenderAxes,'clim',[0 0.45]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
colormap(v.RenderAxes,'gray');
im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
% Partial Effects
for i=1:1:length(Result.PT)
   if isempty(Result.PT(i)), continue; end
    RefScan = clone(RefScan);
    v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
    RefScan.Value = Result.MT.LocalE(4,:,i);set(v.Figure,'Name',[name ' Par ' num2str(i) ' Effect']);set(v.Figure,'Position',get(vref.Figure,'Position'));
    syncCamera(vref,v);
    set(v.RenderAxes,'clim',[0 max(RefScan.Value)]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
    colormap(v.RenderAxes,'gray');
    %set(v.RenderAxes,'clim',[0 35]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
    im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
end
%% PARTIAL R2 & Significance
for T=1:1:length(Result.PT)
   if isempty(Result.PT(T)), continue; end
    RefScan = clone(RefScan);
    v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
    RefScan.Value = Result.PT(T).LocalR2;set(v.Figure,'Name',[name ' Par ' num2str(T) ' R2']);set(v.Figure,'Position',get(vref.Figure,'Position'));
    %set(v.RenderAxes,'clim',[0 max(RefScan.Value)]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
    set(v.RenderAxes,'clim',[0 0.45]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
    colormap(v.RenderAxes,'gray');
    im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
    RefScan = clone(RefScan);
    try
    v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
    RefScan.Value = double((Result.PT(T).pLocalR2<=0.001))';colormap(v.RenderAxes,'summer');set(v.Figure,'Position',get(vref.Figure,'Position'));
    set(v.Figure,'Name',[name ' Par ' num2str(T) ' Sign']);
    colormap(v.RenderAxes,'gray');
    im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
    catch
    end
end

end