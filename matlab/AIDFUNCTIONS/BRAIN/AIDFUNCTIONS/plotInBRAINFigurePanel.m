function fp = plotInBRAINFigurePanel(f,rend,SegmentVal,range,cmap,crit,showcolorbar)
         if nargin<6, crit = []; end
         if nargin<5, cmap = parula(256);end
         if nargin<4, range = [min(SegmentVal), max(SegmentVal)];end
         %if isempty(range), range = [min(SegmentVal), max(SegmentVal)];end
         fp = BRAINFigurePanel(f);
         fp.rosette.select();
         peer = gca;cla(peer);
         plotHierarchicalLabelsDGors(peer,rend,SegmentVal, range, cmap, crit);
         if showcolorbar, colorbar(peer,'Location','SouthOutside');end
         [~,VertexVal] = SegmentVal2VertexVal(rend,SegmentVal);
         for i=1:9
               %i=1;
               fp.ph{i,1}.select();
               peer = gca;cla(peer);
               colormap(peer,cmap);        
               if ~isempty(range), set(peer,'clim',range);end
               renderBrainSurface(rend,VertexVal(i,:),peer);
               view(peer,rend.viewval(1),0);
               light = camlight(peer,'headlight');
               set(light,'Position',get(peer,'CameraPosition'));
               drawnow;
               %title(out.ax1{i},titlenames(levels(i)));
               fp.ph{i,2}.select();
               peer = gca;cla(peer);
               colormap(peer,cmap);        
               if ~isempty(range), set(peer,'clim',range);end
               renderBrainSurface(rend,VertexVal(i,:),peer);
               view(peer,-1*rend.viewval(1),0);
               light = camlight(peer,'headlight');
               set(light,'Position',get(peer,'CameraPosition'));
               drawnow;
         end
end