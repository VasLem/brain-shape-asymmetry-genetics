function out = getValues(obj,uv)
    n = obj.ImageSize(1);
    m = obj.ImageSize(2);
    if max(uv(:))<=1 && min(uv(:))>=0
       %normalised = true;
       x = [0 1];
       y = [0 1];
       %res = [1/n 1/m];
    else
       x = [1 n];
       y = [1 m];
       %res = [1 1];
    end
    xmin = min(x(:)); ymin = min(y(:));
    xmax = max(x(:)); ymax = max(y(:));
    dx = (xmax-xmin)/max(n-1,1);
    dy = (ymax-ymin)/max(m-1,1);
    xx = xmin:dx:xmax;
    yy = ymin:dy:ymax;
    %[xx yy] = meshgrid(0:dx:1,0:dy:1);
    out = zeros(obj.Dim,size(uv,2));
    for i=1:1:obj.Dim
        out(i,:) = interp2(xx,yy,obj.Image(:,:,i),uv(1,:),uv(2,:),obj.PixInterp);
    end
    if obj.Dim==3% RGB values must be within the range of [0 1]
       out = min(ones(size(out)),out);
       out = max(zeros(size(out)),out);
    end
    
end