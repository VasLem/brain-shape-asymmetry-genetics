function [out] = getRIP(in,M)
            out = dot(in',repmat(M/norm(M),1,size(in,1)));
end