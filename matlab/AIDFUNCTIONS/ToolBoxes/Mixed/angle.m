function [out,outd] = angle(v1,v2)
%       if ~size(v1,1)==3, v1 = v1'; end
%       if ~size(v1,1)==3, v2 = v2'; end
%       if ~size(v1,2)==size(v2,2), error('vectors must be same length'); end
%       v1 = v1./norm(v1);
%       v2 = v2./norm(v2);
%       out = acosd(sum(v1.*v2));
      %v1 = v1./norm(v1);
      %v2 = v2./norm(v2);
      T = v1'*v2;
      N = sqrt((v1'*v1)*(v2'*v2));
      out = T/N;
      outd = acosd(out);
end