function imageFacePrediction(obj,RefScanV,vref,name,saveresults)
        BaseFace = obj.BaseFace;
        AmplFace = obj.AmplFace;
        comparison = obj.IdAss;
        BaseFace.Material = 'Dull';
        AmplFace.Material = 'Dull';
%% imaging the scan
        [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Vertices',BaseFace.Vertices,'Texture',[name '_Base'],[]);
        if saveresults,savingResults(v,RefScan);end
%% imaging the norm
        [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Vertices',AmplFace.Vertices,'Texture',[name '_Pred'],[]);
        if saveresults,savingResults(v,RefScan);end
%% imaging the magnitude of the displacements
        val = vDistances(BaseFace,AmplFace);
        [v,RefScan] = cloneAndDisplay(vref,AmplFace,'Value',val,'Indexed',[name 'Disp (Blue NO disp; Red MAX disp)'],[0 max(val)]);
        if saveresults,savingResults(v,RefScan);end
%% imaging the displacements along the normals
        [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Value',comparison.NormalDistances,'Indexed',[name ' Normal (Blue Inward disp; Red Outward disp)'],[]);
        if saveresults,savingResults(v,RefScan);end
%% imaging the Area Ratios
        [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Value',comparison.AreaRatios,'Indexed',[name ' Area (Blue Decrease; Red Increase )'],[]);
        if saveresults,savingResults(v,RefScan);end
%% imaging the signed Curvatures
        maxcurv = min(max(abs(comparison.signedCurvature1)),max(abs(comparison.signedCurvature2)));
        [v,RefScan] = cloneAndDisplay(vref,BaseFace,'Value',comparison.signedCurvature2,'Indexed',[name ' Sign. Curv. Base (Blue Concave; Red Convex)'],[-1*maxcurv maxcurv]);
        if saveresults,savingResults(v,RefScan);end
        [v,RefScan] = cloneAndDisplay(vref,AmplFace,'Value',comparison.signedCurvature1,'Indexed',[name ' Sign. Curv. Pred. (Blue Concave; Red Convex)'],[-1*maxcurv maxcurv]);
        if saveresults,savingResults(v,RefScan);end
%% imaging the signed curvature Differences
        [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Value',comparison.signedCurvatureDiff,'Indexed',[name ' Sign. Curv. Diff (Blue Decreased Convexity; Red Increased convexity)'],[]);
        if saveresults,savingResults(v,RefScan);end
end

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
function savingResults(v,RefScan)
          %save(get(v.Figure,'Name'),'RefScan');
          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
end