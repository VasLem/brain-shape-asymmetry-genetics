function [angle,cosangle] = vectorAngle(v1,v2)
      if ~size(v1,1)==3, v1 = v1'; end
      if ~size(v1,1)==3, v2 = v2'; end
      if ~size(v1,2)==size(v2,2), error('vectors must be same length'); end
      v1 = v1./repmat(sqrt(sum(v1.^2)),3,1);
      v2 = v2./repmat(sqrt(sum(v2.^2)),3,1);
      cosangle = (sum(v1.*v2));
      angle = acosd(cosangle);
end