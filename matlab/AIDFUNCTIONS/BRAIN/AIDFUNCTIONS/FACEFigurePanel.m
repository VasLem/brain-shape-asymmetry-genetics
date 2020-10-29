function fp = FACEFigurePanel(f)
         if nargin<1, f = figure; f.Position = [4 893 1332 433]; end
         clf(f);
         p = panel(f);
         p.pack(1,{1/4 []});
         rosette = p(1,1);
         surface = p(1,2);
         surface.pack(2,3);
         ph = cell(6);
         cols = [1 2 3];
         titlenames = {'i' 'ii' 'iii' 'iv' 'v' 'vi'};
         for i=1:6
             %i=2
             if i<4
                 row = 1;
                 sub = 0;
             elseif i<7
                 row = 2;
                 sub = 3;
             else
                 row = 3;
                 sub = 6;
             end
             colid = i-sub;
             col = cols(colid);
             ph{i} = surface(row,col);
             ph{i}.select();
             title(ph{i},titlenames(i));
         end
         p.select('all');
         p.de.margin = 8;
         rosette.margin = 1;
         fp.f = f;
         fp.p = p;
         fp.rosette = rosette;
         fp.surface = surface;
         fp.ph = ph;
end