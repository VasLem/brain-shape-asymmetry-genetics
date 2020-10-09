function out = getRIP(Ref,std,in,M)
            [~,n] = size(in); % determine input size
            out = nan*zeros(1,n);% allocate memory
            if isempty(std), std = ones(length(Ref),1); end
            % distance between input faces and used reference
            Dir = repmat(Ref,1,n);
            Dir = (in-Dir);
            dist = sqrt(sum(((Dir./repmat(std,1,n)).^2)));
            % Direction between input face and used reference
            coeff2 = M(end,:)'./std;
            for i=1:1:n
                coeff1 = Dir(:,i)/norm(Dir(:,i));
                % Direction of RIP model
                coeff1 = coeff1./std;
                % Angle between direction and model direction
                T = coeff1'*coeff2;
                N = sqrt((coeff1'*coeff1)*(coeff2'*coeff2));
                angle = T/N;
                % parallell distance (RIP)
                out(i) = angle*dist(i);
            end   
end