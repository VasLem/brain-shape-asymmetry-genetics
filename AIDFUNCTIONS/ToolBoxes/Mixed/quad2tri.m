function out = quad2tri(quad)
         Tri1 = quad(1:3,:);
         Tri2 = [quad(3,:);quad(4,:);quad(1,:)];
         out = [Tri1 Tri2];
end