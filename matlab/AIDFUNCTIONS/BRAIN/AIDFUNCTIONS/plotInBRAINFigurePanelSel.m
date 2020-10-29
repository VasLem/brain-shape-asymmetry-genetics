function p = plotInBRAINFigurePanelSel(f,rend,SegmentVal,range,cmap)
         if nargin<5, cmap = parula(256);end
         if nargin<4, range = [min(SegmentVal), max(SegmentVal)];end
         %if isempty(range), range = [min(SegmentVal), max(SegmentVal)];end 
         clf(f);
         p = panel(f);
         p.pack(1,4);
         for i=1:1:4
            pan{i} = p(1,i);
            pan{i}.pack(1,2);
            p1{i} = pan{i}(1,1);
            p2{i} = pan{i}(1,2);
         end
         p.select('all');
         p.margin = 1;
         [~,VertexVal] = SegmentVal2VertexVal(rend,SegmentVal);
         counter = 1;
         for i=2:5
               %i=1;
               p1{counter}.select();
               peer = gca;cla(peer);
               colormap(peer,cmap);        
               if ~isempty(range), set(peer,'clim',range);end
               renderBrainSurface(rend,VertexVal(i,:),peer);
               view(peer,rend.viewval(1),0);
               light = camlight(peer,'headlight');
               set(light,'Position',get(peer,'CameraPosition'));
               drawnow;
               %title(out.ax1{i},titlenames(levels(i)));
               p2{counter}.select();
               peer = gca;cla(peer);
               colormap(peer,cmap);        
               if ~isempty(range), set(peer,'clim',range);end
               renderBrainSurface(rend,VertexVal(i,:),peer);
               view(peer,-1*rend.viewval(1),0);
               light = camlight(peer,'headlight');
               set(light,'Position',get(peer,'CameraPosition'));
               drawnow;
               counter = counter+1;
         end
end