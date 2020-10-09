function fp = BRAINFigurePanel(f)
         if nargin<1, f = figure; f.Position = [4 893 1332 433]; end
         clf(f);
         p = panel(f);
         p.pack(1,{1/4 []});
         rosette = p(1,1);
         surface = p(1,2);
         surface.pack(3,6);
         ph = cell(9,2);
         cols = [1 2;3 4;5 6];
         titlenames = {'i' 'ii' 'iii' 'iv' 'v' 'vi' 'vii' 'viii' 'ix'};
         for i=1:9
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
             col = cols(colid,:);
             ph{i,1} = surface(row,col(1));
             ph{i,1}.select();
             title(ph{i,1},titlenames(i));
             ph{i,2} = surface(row,col(2));
         end
         p.select('all');
         p.de.margin = 4;
         rosette.margin = 1;
         fp.f = f;
         fp.p = p;
         fp.rosette = rosette;
         fp.surface = surface;
         fp.ph = ph;
end