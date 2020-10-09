function out = readXMLLandmarks(file)
         A=fileread(file);
         B=strfind(A,'<pos>'); 
         E=strfind(A,'</pos>');
         n = length(B);
         Vertices = zeros(n,3);
         for i=1:1:n
             %i=1;
             tmp = A(B(i)+5:E(i)-1);
             index = strfind(tmp,' ');
             Vertices(i,1) = str2double(tmp(1:index(1)-1));
             Vertices(i,2) = str2double(tmp(index(1)+1:index(2)-1));
             Vertices(i,3) = str2double(tmp(index(2)+1:index(3)-1));
         end
         out = shape3D;
         out.Vertices = Vertices;
end