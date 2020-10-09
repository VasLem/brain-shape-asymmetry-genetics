function showAnovaResults(A, scan, prefix, vref)
         scan = clone(scan);
         scan.Value = ones(1,scan.nrV);
         scan.ColorMode = 'Indexed';
         scan.Material = 'Dull';
         if nargin < 4
            vref = viewer(scan);
            pause;
            disp('press enter to continue');
         end
         % Individual 
             showValResult(A.LM.I,scan,vref,[prefix ' Ind MS ']);%num2str(A.Total.I)]);
             showValResult(A.LM.IF,scan,vref,[prefix ' Ind F ']);%num2str(A.Total.IF)]);
             try
                showPResult(A.LM.permIF,scan,vref,[prefix ' Ind P ']);% num2str(A.Total.IP)]);
             catch
                showPResult(A.LM.IP,scan,vref,[prefix ' Ind P ']);% num2str(A.Total.IP)]);
             end
         % Directional 
             showValResult(A.LM.D,scan,vref,[prefix ' Dir MS ']);% num2str(A.Total.D)]);
             showValResult(A.LM.DF,scan,vref,[prefix ' Dir F ']);% num2str(A.Total.DF)]);
             try
                showPResult(A.LM.permDF,scan,vref,[prefix ' Dir P ']);% num2str(A.Total.DP)]);
             catch
                showPResult(A.LM.DP,scan,vref,[prefix ' Dir P ']);% num2str(A.Total.DP)]);
             end
         % Fluctuating 
             showValResult(A.LM.F,scan,vref,[prefix ' Fluc MS ']);% num2str(A.Total.F)]);
             showValResult(A.LM.FF,scan,vref,[prefix ' Fluc F ']);% num2str(A.Total.FF)]);
             try
                showPResult(A.LM.permFF,scan,vref,[prefix ' Fluc P ']);% num2str(A.Total.FP)]);
             catch
                showPResult(A.LM.FP,scan,vref,[prefix ' Fluc P ']);% num2str(A.Total.FP)]);
             end
         % Error
             showValResult(A.LM.E,scan,vref,[prefix ' Err MS ']);% num2str(A.Total.E)]);
end

function showValResult(val,scan,vref,name)
         s = clone(scan);
         s.Value = val;
         v = viewer(s);
         if ~isempty(vref), syncCamera(vref,v); end
         v.BackgroundColor = [1 1 1];
         set(v.Figure,'Name',name);
         colormap(v.RenderAxes,'jet');
         colorbar('peer',v.RenderAxes);
         set(v.Toolbar.light_toggle,'State','on');
         set(v.Toolbar.link_toggle,'State','on');
         if ~isempty(vref), copyViewer(vref,v); end
end

function showPResult(val,scan,vref,name)
         tmp = ones(size(val));
         tmp(val>=0.001) = 0.5;
         tmp(val>=0.05) = 0;
         s = clone(scan);
         s.Value = tmp;
         s.Material = 'Facial';
         v = viewer(s);
         syncCamera(vref,v);
         v.BackgroundColor = [1 1 1];
         set(v.RenderAxes,'clim',[0 1]);
         colormap(v.RenderAxes,'summer');
         set(v.Figure,'Name',name);
         set(v.Toolbar.light_toggle,'State','on');
         set(v.Toolbar.link_toggle,'State','on');
end