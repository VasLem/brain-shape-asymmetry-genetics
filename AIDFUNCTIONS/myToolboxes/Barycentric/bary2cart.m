function [cart] = bary2cart(Vertices,Faces,bar,Findex)
        cart = zeros(size(bar));
        for i=1:1:length(Findex)
            f = Findex(i);
            fp = Vertices(:,Faces(:,f));
            a = fp(:,1);b=fp(:,2);c=fp(:,3);
            cart(:,i) = bar(1,i)*a + bar(2,i)*b + bar(3,i)*c;
        end
end