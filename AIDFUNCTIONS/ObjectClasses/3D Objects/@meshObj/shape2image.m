function out = shape2image(obj,res,range)
        if nargin<3
           %range(:,1) = min(obj.Vertices');
           %range(:,2) = max(obj.Vertices');
        end
        UV = obj.UV;
        for i=1:1:2
            UV(i,:) = UV(i,:)-min(UV(i,:));
            UV(i,:) = UV(i,:)/max(UV(i,:));
        end
        %res = [0.01 0.01];
        %if isscalar(res), res = [res res]; end
        res = (1/(res-1));
        [X,Y] = meshgrid(0:res:1,0:res:1);
        im = zeros(size(X,1),size(X,2),3);
        warning('off','All');
        for i=1:1:3
            im(:,:,i) = griddata(UV(1,:),UV(2,:),obj.Vertices(i,:),X,Y,'cubic');
            im(:,:,i) = (im(:,:,i)-range(i,1))/(range(i,2)-range(i,1));
        end
        im = im*255;
        out = uint8(im);
end