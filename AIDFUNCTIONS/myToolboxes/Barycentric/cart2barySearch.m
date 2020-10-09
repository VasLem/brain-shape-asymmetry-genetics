function [bar,index] = cart2barySearch(MeshVertices,cart)
        index = knnsearch(MeshVertices,cart,'K',3);% look for three closest points
        bar = zeros(size(cart));
        n = size(cart,1);
        for i=1:1:length(Findex)
            f = Findex(i);
            p = cart(:,i);
            fp = Vertices(:,Faces(:,f));
            a = fp(:,1);b=fp(:,2);c=fp(:,3);
            n = cross((b-a),(c-a));
            na = cross((c-b),(p-b));
            nb = cross((a-c),(p-c));
            nc = cross((b-a),(p-a));
            bar(1,i) = (dot(n,na))/norm(n)^2;
            bar(2,i) = (dot(n,nb))/norm(n)^2;
            bar(3,i) = (dot(n,nc))/norm(n)^2;
        end
end