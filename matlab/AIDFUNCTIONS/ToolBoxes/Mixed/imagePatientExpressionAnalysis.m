function imagePatientExpressionAnalysis(res,RefScanV,vref,name,saveresults)
res.assorig.Norm.Material = 'Dull';
delete(res.assorig.Norm.Border);
delete(res.assorig.Norm.PoseLM);
res.assorig.Scan.Material = 'Dull';
delete(res.assorig.Scan.Border);
delete(res.assorig.Scan.PoseLM);
res.assrefl.Scan.Material = 'Dull';
delete(res.assrefl.Scan.Border);
delete(res.assrefl.Scan.PoseLM);
[v,RefScan] = cloneAndDisplay(vref,RefScanV,'Vertices',res.assorig.Scan.Vertices,'Single',[name ' repose'],[]);
if saveresults,savingResults(v,RefScan);end
[v,RefScan] = cloneAndDisplay(vref,RefScanV,'Vertices',res.assorig.Norm.Vertices,'Single',[name ' smile'],[]);
if saveresults,savingResults(v,RefScan);end
[v,RefScan] = cloneAndDisplay(vref,RefScanV,'Vertices',res.assrefl.Scan.Vertices,'Single',[name ' repose Refl.'],[]);
if saveresults,savingResults(v,RefScan);end
[v,RefScan] = cloneAndDisplay(vref,RefScanV,'Vertices',res.assrefl.Norm.Vertices,'Single',[name ' smile Refl.'],[]);
if saveresults,savingResults(v,RefScan);end
%% imaging the expression distanceMap
[v,RefScan] = cloneAndDisplay(vref,res.assorig.Norm,'Value',res.assorig.DistanceMap,'Indexed',[name 'Dist. Expr.'],[0 10]);
if saveresults,savingResults(v,RefScan);end
%% imaging the reflected expression distanceMap
res.assrefl.Norm.Material = 'Dull';
delete(res.assrefl.Norm.Border);
delete(res.assrefl.Norm.PoseLM);
[v,RefScan] = cloneAndDisplay(vref,res.assrefl.Norm,'Value',res.assrefl.DistanceMap,'Indexed',[name 'Dist. Expr. Refl.'],[0 10]);
if saveresults,savingResults(v,RefScan);end
%% image the asymmetry in the distancemap
res.SymSmilePatient.Material = 'Dull';
delete(res.SymSmilePatient.Border);
delete(res.SymSmilePatient.PoseLM);
[v,RefScan] = cloneAndDisplay(vref,res.SymSmilePatient,'Value',res.AsymmetryDistance,'Indexed',[name 'Dist. Expr. Asym.'],[0 10]);
if saveresults,savingResults(v,RefScan);end
%% image the curvature onto the repose image
val = max(nanmax(abs([res.Curvature.cneutral res.Curvature.csmile])));
[v,RefScan] = cloneAndDisplay(vref,res.assorig.Scan,'Value',res.Curvature.cneutral,'Indexed',[name 'Curv. Repose'],[-val val]);
if saveresults,savingResults(v,RefScan);end
[v,RefScan] = cloneAndDisplay(vref,res.assorig.Norm,'Value',res.Curvature.csmile,'Indexed',[name 'Curv. Smile'],[-val val]);
if saveresults,savingResults(v,RefScan);end
%% image the curvature change due to the expression
[v,RefScan] = cloneAndDisplay(vref,res.assorig.Norm,'Value',res.Curvature.difforig,'Indexed',[name 'Curv. Expr.'],[0 0.15]);
if saveresults,savingResults(v,RefScan);end
%% image the curvature change due to the expression
[v,RefScan] = cloneAndDisplay(vref,res.assrefl.Norm,'Value',res.Curvature.diffrefl,'Indexed',[name 'Curv. Expr. Refl.'],[0 0.15]);
if saveresults,savingResults(v,RefScan);end
%% image the area change due to the expression
[v,RefScan] = cloneAndDisplay(vref,res.assorig.Norm,'Value',res.Area.origRatios,'Indexed',[name 'Area Expr.'],[-1 1]);
if saveresults,savingResults(v,RefScan);end
%% image the area change due to the expression
[v,RefScan] = cloneAndDisplay(vref,res.assrefl.Norm,'Value',res.Area.reflRatios,'Indexed',[name 'Area Expr. Refl.'],[-1 1]);
if saveresults,savingResults(v,RefScan);end
%% image the normal displacement due to the expression
[v,RefScan] = cloneAndDisplay(vref,res.assorig.Norm,'Value',res.Normal.ndist,'Indexed',[name 'Norm Disp. Expr.'],[-10 10]);
if saveresults,savingResults(v,RefScan);end
%% image the normal displacement due to the expression
[v,RefScan] = cloneAndDisplay(vref,res.assrefl.Norm,'Value',res.Normal.ndistRefl,'Indexed',[name 'Norm Disp. Expr. Refl.'],[-10 10]);
if saveresults,savingResults(v,RefScan);end
%% image the normal displacement due to the expression
[v,RefScan] = cloneAndDisplay(vref,res.SymSmilePatient,'Value',res.Normal.DIFF,'Indexed',[name 'Norm Disp. Expr. Asym.'],[-10 10]);
if saveresults,savingResults(v,RefScan);end
end


function [v,RefScan] = cloneAndDisplay(vref,RefScanV,what,value,color,name,range)
         if nargin < 7, range = []; end
         RefScan = clone(RefScanV);delete(RefScan.PoseLM);
         v = viewer(RefScan);v.BackgroundColor = [0 0 0];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
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
         end
         set(v.Figure,'Name',name);set(v.Figure,'Position',get(vref.Figure,'Position'));
         v.SceneLightPosition = vref.SceneLightPosition;
         drawnow;
end
function savingResults(v,RefScan)
          %save(get(v.Figure,'Name'),'RefScan');
          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
end