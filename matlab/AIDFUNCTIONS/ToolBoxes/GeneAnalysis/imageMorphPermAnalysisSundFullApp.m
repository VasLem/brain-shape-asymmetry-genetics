function imageMorphPermAnalysisSundFullApp(analysis,RefScanV,vref,name,mi,pl,saveresults,gcode,texturescan)
%% imaging the scan
        [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Vertices',analysis.Morph2.Vertices,'Single',[name ' ' pl],[]);
        if ~isempty(analysis.Morph2.PoseLM)
           RefScan.PoseLM = clone(analysis.Morph2.PoseLM);
           RefScan.PoseLM.Visisble = true;
        end
        switch gcode
            case -1
                v.BackgroundColor = [0.95 0.87 0.73];
            case 1
                v.BackgroundColor = [0.76 0.87 0.78];
            otherwise
                v.BackgroundColor = [1 1 1];
        end
        RefScan.UV = texturescan.UV;
        RefScan.TextureMap = clone(texturescan.TextureMap);
        RefScan.ColorMode = 'Texture';
        if saveresults,savingResults(v,RefScan);end
        morph2 = clone(RefScan);
%% imaging the norm
        [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Vertices',analysis.Morph1.Vertices,'Single',[name ' ' mi],[]);
        if ~isempty(analysis.Morph1.PoseLM)
           RefScan.PoseLM = clone(analysis.Morph1.PoseLM);
           RefScan.PoseLM.Visisble = true;
        end
        switch gcode
            case 1
                v.BackgroundColor = [0.95 0.87 0.73];
            case -1
                v.BackgroundColor = [0.76 0.87 0.78];
            otherwise
                v.BackgroundColor = [1 1 1];
        end
        RefScan.UV = texturescan.UV;
        RefScan.TextureMap = clone(texturescan.TextureMap);
        RefScan.ColorMode = 'Texture';
        if saveresults,savingResults(v,RefScan);end
        morph1 = clone(RefScan);
%% imaging the effect map
        [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Value',analysis.STATS.LocalE(4,:),'Indexed',[name ' Effect'],[0 max(analysis.STATS.LocalE(4,:))]);
        if saveresults,savingResults(v,RefScan);end
% %% imaging the effect-size map
%         [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Value',analysis.STATS.LocalS,'Indexed',[name ' Effect Size'],[0 max(analysis.STATS.LocalS)]);
%         if saveresults,savingResults(v,RefScan);end
% %% imaging the significance of the effect
%         [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Value',double(analysis.STATS.pLocalS<=0.001)','Indexed',[name ' Effect Sign.'],[0 1]);
%         colormap(v.RenderAxes,'summer');colorbar('peer',v.RenderAxes,'off');
%         if saveresults,savingResults(v,RefScan);end
% %% imaging the X disp
%         [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Value',analysis.STATS.comparison.Differences(1,:),'Indexed',[name ' X (Blue Left disp; Red Right disp)'],[]);
%         if saveresults,savingResults(v,RefScan);end
% %% imaging the Y disp
%         [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Value',analysis.STATS.comparison.Differences(2,:),'Indexed',[name 'Y (Blue Downward disp; Red Upward disp)'],[]);
%         if saveresults,savingResults(v,RefScan);end   
% %% imaging the Z disp
%         [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Value',analysis.STATS.comparison.Differences(3,:),'Indexed',[name 'Z (Blue Backwards disp; Red Forward disp)'],[]);
%         if saveresults,savingResults(v,RefScan);end        
%% imaging the displacements along the normals
        [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Value',analysis.STATS.comparison.NormalDistances,'Indexed',[name ' Normal (Blue Inward disp; Red Outward disp)'],[]);
        if saveresults,savingResults(v,RefScan);end
% %% imaging H1 pos of the displacements along the normals
%         [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Value',double(analysis.STATS.comparison.pH1PosNormalDistances<=0.001)','Indexed',[name ' Normal H1 +'],[0 1]);
%         colormap(v.RenderAxes,'summer');colorbar('peer',v.RenderAxes,'off');
%         if saveresults,savingResults(v,RefScan);end
% %% imaging H1 neg of the displacements along the normals
%         [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Value',double(analysis.STATS.comparison.pH1NegNormalDistances<=0.001)','Indexed',[name ' Normal H1 -'],[0 1]);
%         colormap(v.RenderAxes,'summer');colorbar('peer',v.RenderAxes,'off');
%         if saveresults,savingResults(v,RefScan);end 
% %% imaging H2 of the displacements along the normals
%         [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Value',double(analysis.STATS.comparison.pH2NormalDistances<=0.001)','Indexed',[name ' Normal H2'],[0 1]);
%         colormap(v.RenderAxes,'summer');colorbar('peer',v.RenderAxes,'off');
%         if saveresults,savingResults(v,RefScan);end
%% imaging the Area Ratios
        [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Value',analysis.STATS.comparison.AreaRatios,'Indexed',[name ' Area (Blue Decrease; Red Increase )'],[]);
        if saveresults,savingResults(v,RefScan);end
% %% imaging H1 pos Area Ratios
%         [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Value',double(analysis.STATS.comparison.pH1PosAreaRatios<=0.001)','Indexed',[name ' Area H1 +'],[0 1]);
%         colormap(v.RenderAxes,'summer');colorbar('peer',v.RenderAxes,'off');
%         if saveresults,savingResults(v,RefScan);end
% %% imaging H1 neg Area Ratios
%         [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Value',double(analysis.STATS.comparison.pH1NegAreaRatios<=0.001)','Indexed',[name ' Area H1 -'],[0 1]);
%         colormap(v.RenderAxes,'summer');colorbar('peer',v.RenderAxes,'off');
%         if saveresults,savingResults(v,RefScan);end 
% %% imaging H2 Area Ratios
%         [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Value',double(analysis.STATS.comparison.pH2AreaRatios<=0.001)','Indexed',[name ' Area H2'],[0 1]);
%         colormap(v.RenderAxes,'summer');colorbar('peer',v.RenderAxes,'off');
%         if saveresults,savingResults(v,RefScan);end
%% imaging the signed Curvatures
        maxcurv = min(max(abs(analysis.STATS.comparison.signedCurvature1)),max(abs(analysis.STATS.comparison.signedCurvature2)));
        [v,RefScan] = cloneAndDisplay(vref,morph2,'Value',analysis.STATS.comparison.signedCurvature2,'Indexed',[name ' Sign. Curv. ' pl ' (Blue Concave; Red Convex)'],[-1*maxcurv maxcurv]);
        if saveresults,savingResults(v,RefScan);end
        [v,RefScan] = cloneAndDisplay(vref,morph1,'Value',analysis.STATS.comparison.signedCurvature1,'Indexed',[name ' Sign. Curv. ' mi ' (Blue Concave; Red Convex)'],[-1*maxcurv maxcurv]);
        if saveresults,savingResults(v,RefScan);end
% %% imaging the signed curvature Differences
%         [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Value',analysis.STATS.comparison.signedCurvatureDiff,'Indexed',[name ' Sign. Curv. Diff (Blue Decreased Convexity; Red Increased convexity)'],[]);
%         if saveresults,savingResults(v,RefScan);end
% %% imaging H1 pos signed curvature differences
%         [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Value',double(analysis.STATS.comparison.pH1PosCurvDiff<=0.001)','Indexed',[name ' Sign. Curv. H1 +'],[0 1]);
%         colormap(v.RenderAxes,'summer');colorbar('peer',v.RenderAxes,'off');
%         if saveresults,savingResults(v,RefScan);end
% %% imaging H1 neg signed curvature differences
%         [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Value',double(analysis.STATS.comparison.pH1NegCurvDiff<=0.001)','Indexed',[name ' Sign. Curv. H1 -'],[0 1]);
%         colormap(v.RenderAxes,'summer');colorbar('peer',v.RenderAxes,'off');
%         if saveresults,savingResults(v,RefScan);end 
% %% imaging H2 signed curvature differences
%         [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Value',double(analysis.STATS.comparison.pH2CurvDiff<=0.001)','Indexed',[name ' Sign. Curv. H2'],[0 1]);
%         colormap(v.RenderAxes,'summer');colorbar('peer',v.RenderAxes,'off');
%         if saveresults,savingResults(v,RefScan);end             
%% imaging the curvature ratios, I don't beleive they give more insight compared to signed curvature
%         [v,RefScan] = cloneAndDisplay(vref,RefScanV,'Value',analysis.STATS.comparison.CurvatureRatios,'Indexed',[name ' Curv. (Blue Decrease; Red Increase)'],[]);
%         if saveresults,savingResults(v,RefScan);end

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
         end
         set(v.Figure,'Name',name);set(v.Figure,'Position',get(vref.Figure,'Position'));
         v.SceneLightPosition = vref.SceneLightPosition;
         drawnow;
end
function savingResults(v,RefScan)
          save(get(v.Figure,'Name'),'RefScan');
          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
end