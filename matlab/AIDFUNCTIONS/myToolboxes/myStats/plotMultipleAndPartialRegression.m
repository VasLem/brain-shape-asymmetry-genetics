function plotMultipleAndPartialRegression(Result,RefScan)
RefScan.Material = 'Facial';
RefScan.Value = ones(1,RefScan.nrV);
RefScan.SingleColor = [0.8 0.8 0.8];
RefScan.ColorMode = 'Indexed';
vref = viewer(RefScan);
vref.BackgroundColor = [1 1 1];
vref.SceneLightLinked = true;
vref.SceneLightVisible = true;
pause;
%% Multiple effect
RefScan = clone(RefScan);
v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
RefScan.Value = double((Result.MT.pLocalR2<=0.001))';colormap(v.RenderAxes,'summer');
set(v.Figure,'Name','Mult. Sign.');
RefScan = clone(RefScan);
v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
RefScan.Value = Result.MT.LocalR2;set(v.Figure,'Name','Mult. R2');
%% Partial Effects
for i=1:1:length(Result.PT)
    if isempty(Result.PT(i)), continue; end
    RefScan = clone(RefScan);
    v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
    RefScan.Value = Result.MT.LocalE(4,:,i);set(v.Figure,'Name',['Par. ' num2str(i) ' Effect']);
    syncCamera(vref,v);
end
%% PARTIAL R2 & Significance
for T=1:1:length(Result.PT)
    if isempty(Result.PT(T)), continue; end
    RefScan = clone(RefScan);
    v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
    RefScan.Value = Result.PT(T).LocalR2;set(v.Figure,'Name',['Par. ' num2str(T) ' R2']);
    RefScan = clone(RefScan);
    v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
    RefScan.Value = double((Result.PT(T).pLocalR2<=0.001))';colormap(v.RenderAxes,'summer');
    set(v.Figure,'Name',['Par. ' num2str(T) ' Sign.']);
end

end