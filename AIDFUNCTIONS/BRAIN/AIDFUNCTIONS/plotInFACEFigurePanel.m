function fp = plotInFACEFigurePanel(f,rend,SegmentVal,range,cmap,crit,showcolorbar)
         if nargin<6, crit = []; end
         if nargin<5, cmap = parula(256);end
         if nargin<4, range = [min(SegmentVal), max(SegmentVal)];end
         %if isempty(range), range = [min(SegmentVal), max(SegmentVal)];end
         fp = FACEFigurePanel(f);
         fp.rosette.select();
         peer = gca;cla(peer);
         plotHierarchicalLabelsDGors(peer,rend,SegmentVal, range, cmap, crit);
         if showcolorbar, colorbar(peer,'Location','SouthOutside');end
         [~,VertexVal] = SegmentVal2VertexVal(rend,SegmentVal);
         for i=1:6
               %i=1;
               fp.ph{i}.select();
               peer = gca;cla(peer);
               colormap(peer,cmap);        
               if ~isempty(range), set(peer,'clim',range);end
               renderFaceSurface(rend,VertexVal(i,:),peer);
               light = camlight(peer,'headlight');
               set(light,'Position',get(peer,'CameraPosition'));
               drawnow;
         end
end