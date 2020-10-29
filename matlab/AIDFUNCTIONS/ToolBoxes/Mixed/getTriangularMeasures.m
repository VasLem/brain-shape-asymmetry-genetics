function [measures, measurenames] = getTriangularMeasures(in)
         P1 = in(:,1);
         P2 = in(:,2);
         P3 = in(:,3);
         D1 = sqrt(sum((P1-P2).^2));
         D2 = sqrt(sum((P1-P3).^2));
         D3 = sqrt(sum((P2-P3).^2));
         A1 = acosd((-D3^2+D1^2+D2^2)/(2*D1*D2));
         A2 = acosd((D2^2-D1^2-D3^2)/(-2*D1*D3));
         A3 = acosd((D1^2-D2^2-D3^2)/(-2*D2*D3));
         measures = [P1;P2;P3;D1;D2;D3;A1;A2;A3];
         measurenames = {'P1X' 'P1Y' 'P1Z' 'P2X' 'P2Y' 'P2Z'...
                         'P3X' 'P3Y' 'P3Z' 'D1' 'D2' 'D3' 'A1' 'A2' 'A3'};
end