function imageMorphAssessment(ass,RefScanV,vref,name)
         % imaging the scan
         RefScan = clone(RefScanV);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Vertices = ass.Norm.Vertices;
         RefScan.ColorMode = 'Single';
         set(v.Figure,'Name',[name ' AA']);set(v.Figure,'Position',get(vref.Figure,'Position'));
         v.SceneLightPosition = vref.SceneLightPosition;
         drawnow;
%          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
         % imaging the norm
         RefScan = clone(RefScanV);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Vertices = ass.Scan.Vertices;
         RefScan.ColorMode = 'Single';
         set(v.Figure,'Name',[name ' BB']);set(v.Figure,'Position',get(vref.Figure,'Position'));
         v.SceneLightPosition = vref.SceneLightPosition;
         drawnow;
%          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
         % imaging the distance map
         RefScan = clone(RefScanV);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = ass.DistanceMap;
         RefScan.ColorMode = 'Indexed';
         set(v.RenderAxes,'clim',[0 max(RefScan.Value)]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name ' DistanceMap']);set(v.Figure,'Position',get(vref.Figure,'Position'));
%          im = captureImage(v);imwrite(im,[get(v.Figure,'Name') '.tiff'],'tiff','Compression','none','Resolution',600);drawnow;delete(v);
%          % X disp
         diff = ass.Scan.Vertices-ass.Norm.Vertices;
         RefScan = clone(RefScanV);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = diff(1,:);
         RefScan.ColorMode = 'Indexed';
         set(v.RenderAxes,'clim',[-1*max(abs(RefScan.Value)) max(abs(RefScan.Value))]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name ' X (Blue Left disp; Red Right disp)']);set(v.Figure,'Position',get(vref.Figure,'Position'));
%          % Y disp
         RefScan = clone(RefScanV);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = diff(2,:);
         RefScan.ColorMode = 'Indexed';
         set(v.RenderAxes,'clim',[-1*max(abs(RefScan.Value)) max(abs(RefScan.Value))]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name 'Y (Blue Downward disp; Red Upward disp)']);set(v.Figure,'Position',get(vref.Figure,'Position'));
%          % Z disp
         RefScan = clone(RefScanV);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = diff(3,:);
         RefScan.ColorMode = 'Indexed';
         set(v.RenderAxes,'clim',[-1*max(abs(RefScan.Value)) max(abs(RefScan.Value))]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name 'Z (Blue Backwards disp; Red Forward disp)']);set(v.Figure,'Position',get(vref.Figure,'Position'));
%          % Normal Disp.
         VF = vectorField3D;
         VF.StartPoints = ass.Norm.Vertices;
         VF.EndPoints = ass.Scan.Vertices;
         [~,cosangle] = vectorAngle(-1*ass.Norm.Gradient,VF.Direction);
         RefScan = clone(RefScanV);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = VF.Length.*cosangle;
         RefScan.ColorMode = 'Indexed';
         set(v.RenderAxes,'clim',[-1*max(abs(RefScan.Value)) max(abs(RefScan.Value))]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name ' Normal (Blue Inward disp; Red Outward disp)']);set(v.Figure,'Position',get(vref.Figure,'Position'));
%% Area Changes
          obj = ass.Norm;
          face = obj.Faces;
          nrF = size(face,2);
          LOC =  zeros(3,nrF,3);AB = zeros(3,nrF);AC = zeros(3,nrF);
          for i=1:1:3
                   LOC(:,:,i) = reshape(obj.Vertices(i,face(:)),3,size(face,2));
                   AB(i,:) = LOC(1,:,i)-LOC(2,:,i);
                   AC(i,:) = LOC(1,:,i)-LOC(3,:,i);
          end
          areasNorm = 0.5*sqrt(dot(AB,AB).*dot(AC,AC)-dot(AB,AC).^2);
          avgareasNorm = zeros(1,obj.nrV);
          parfor i=1:obj.nrV
                 tmp = ismember(obj.Faces,i);
                 [~,j] = find(tmp==1);
                 avgareasNorm(i) = mean(areasNorm(j));
          end
          obj = ass.Scan;
          face = obj.Faces;
          nrF = size(face,2);
          LOC =  zeros(3,nrF,3);AB = zeros(3,nrF);AC = zeros(3,nrF);
          for i=1:1:3
                   LOC(:,:,i) = reshape(obj.Vertices(i,face(:)),3,size(face,2));
                   AB(i,:) = LOC(1,:,i)-LOC(2,:,i);
                   AC(i,:) = LOC(1,:,i)-LOC(3,:,i);
          end
          areasScan = 0.5*sqrt(abs(dot(AB,AB).*dot(AC,AC)-dot(AB,AC).^2));
          avgareasScan = zeros(1,obj.nrV);
          parfor i=1:obj.nrV
                 tmp = ismember(obj.Faces,i);
                 [~,j] = find(tmp==1);
                 avgareasScan(i) = mean(areasScan(j));
          end
         
          
         ratios = avgareasNorm./avgareasScan;
         %ratios = -1*(ratios-1);
          
         RefScan = clone(RefScanV);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = -1*log(ratios);
         RefScan.ColorMode = 'Indexed';
         set(v.RenderAxes,'clim',[-1*max(abs(RefScan.Value)) max(abs(RefScan.Value))]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name ' Areal (Blue Decrease; Red Increase )']);set(v.Figure,'Position',get(vref.Figure,'Position'));
%% Curvature Changes
         [Cmean,Cgaussian]=curvature(ass.Norm,true);
         [Cmean2,Cgaussian2]=curvature(ass.Scan,true);
         
         ratiosM = Cmean./Cmean2;
         
         RefScan = clone(RefScanV);
         v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
         RefScan.Value = -1*log(ratiosM);
         RefScan.ColorMode = 'Indexed';
         set(v.RenderAxes,'clim',[-1*max(abs(RefScan.Value)) max(abs(RefScan.Value))]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
         set(v.Figure,'Name',[name ' Mean Curv. (Blue Decrease; Red Increase )']);set(v.Figure,'Position',get(vref.Figure,'Position'));
         
%          ratiosG = Cgaussian./Cgaussian2;
%          
%          RefScan = clone(RefScanV);
%          v = viewer(RefScan);v.BackgroundColor = [1 1 1];v.SceneLightLinked = true;v.SceneLightVisible = true;syncCamera(vref,v);
%          RefScan.Value = -1*log(ratiosG);
%          RefScan.ColorMode = 'Indexed';
%          set(v.RenderAxes,'clim',[-1*max(abs(RefScan.Value)) max(abs(RefScan.Value))]);colorbar('peer',v.RenderAxes,'LOCATION','EastOutside');
%          set(v.Figure,'Name',[name ' Gaussian Curv. (Blue Decrease; Red Increase )']);set(v.Figure,'Position',get(vref.Figure,'Position'));
%          
         
         
    
          
          
         
         
end