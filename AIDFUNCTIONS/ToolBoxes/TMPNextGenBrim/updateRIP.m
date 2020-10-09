function [out] = updateRIP(in,M)
            % in this implementation I take the reference as the origin
            n2 = size(in,1); % determine input size
            in = in';
            n1 = size(M,1);
            out = nan*zeros(n1,n2);% allocate memory
            for j=1:n1
                coeff2 = M(j,:)';
                coeff2 = coeff2/norm(coeff2);
                out(j,:) = dot(in,repmat(coeff2,1,n2));
            end
end