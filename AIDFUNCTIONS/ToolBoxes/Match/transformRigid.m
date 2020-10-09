function points = transformRigid(points,T)
         if ~size(T,1)==4||~size(T,2)==4, error('Incorrect Transformation Matrix'); end
         if ~size(points,1)==3
             points = points';
         end
         points = [points;ones(1,size(points,2))];
         points = T*points;
         points = points(1:3,:);
end