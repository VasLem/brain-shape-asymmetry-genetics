function screePlot(obj, percvar)
         if nargin < 2
            percvar = 90;
         end
         figure;title('Scree Plot');
         X = (1:obj.nrEV);
         Y1 = obj.Explained;
         Y2 = zeros(1,obj.nrEV);
         Y2(1) = Y1(1);
         for i=2:1:obj.nrEV;
             Y2(i) = Y2(i-1)+Y1(i);
         end
         [ax,h1,h2] = plotyy(X,Y1,X,Y2);
         set(ax,'ylim',[0 100]);
         set(ax,'xlim',[X(1) X(end)]);
end