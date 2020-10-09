function points = transformRBF(points,T)
         if ~(size(points,1)==3||size(points,1)==2)
             points = points';
         end
         dim = size(points,1);
         pl = fastrbf_makepointlist(points);clear points;
         plx = fastrbf_pointeval(T.X,pl,'messages',0);
         ply = fastrbf_pointeval(T.Y,pl,'messages',0);
         points = [plx.Value;ply.Value];
         clear plx ply;
         if dim == 3
            plz = fastrbf_pointeval(T.Z,pl,'messages',0);
            points = [points;plz.Value];
            clear plz;
         end
         clear pl;
end