function [cart] = bary2cartKNN(MeshVertices,index,bar)
        bar = bar';
        cart = zeros(size(bar));
        n = size(cart,2);
        for i=1:1:n
            a = MeshVertices(index(i,1),:)';
            b = MeshVertices(index(i,2),:)';
            c = MeshVertices(index(i,3),:)';
%             f = Findex(i);
%             fp = Vertices(:,Faces(:,f));
%             a = fp(:,1);b=fp(:,2);c=fp(:,3);
            cart(:,i) = bar(1,i)*a + bar(2,i)*b + bar(3,i)*c;
        end
        cart = cart';
end