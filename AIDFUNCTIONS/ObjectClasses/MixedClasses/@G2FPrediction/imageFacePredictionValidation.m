function imageFacePredictionValidation(obj,RefScanV,truth,car,opp,vref,name,saveresults)
        BaseFace = obj.BaseFace;
        PredFace = obj.PredFace;
        AmplFace = obj.AmplFace;
        comparison = obj.IdAss;
        BaseFace.Material = 'Dull';
        AmplFace.Material = 'Dull';
%% imaging the truth
        [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Vertices',truth.Vertices,'Single',[name '_truth'],[]);
        if saveresults,savingResults(v,RefScan);end
%% imaging the Car
        [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Vertices',car.Vertices,'Single',[name '_Car'],[]);
        if saveresults,savingResults(v,RefScan);end
%% imaging the Opp
        [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Vertices',opp.Vertices,'Single',[name '_Opp'],[]);
        if saveresults,savingResults(v,RefScan);end        
%% imaging the Base Face
        [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Vertices',BaseFace.Vertices,'Single',[name '_Base'],[]);
        if saveresults,savingResults(v,RefScan);end
%% imaging the Amplified face
        [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Vertices',AmplFace.Vertices,'Single',[name '_Ampl'],[]);
        if saveresults,savingResults(v,RefScan);end
%% imaging the norm
        [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Vertices',PredFace.Vertices,'Single',[name '_Pred'],[]);
        if saveresults,savingResults(v,RefScan);end        
%% imaging the displacements along the normals
        [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Value',comparison.NormalDistances,'Indexed',[name ' Normal (Blue Inward disp; Red Outward disp)'],[]);
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