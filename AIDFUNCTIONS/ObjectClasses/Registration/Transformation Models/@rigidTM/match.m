function out = match(obj,q,p,w)
% find the rigid transform from p to q
         if nargout == 1
            obj = clone(obj);
            out = obj;
         end
         p = TM.getPoints(p);
         q = TM.getPoints(q);
         nrpoints = size(p,2);
         if nargin<4
            w  = ones(3,nrpoints);
         else
            if size(w,1)==1, w = repmat(w,3,1);end 
         end
         totalw = sum(w,2);
         index = find(sum(w));% find points not having weights == 0;
         if length(index)<nrpoints
             p = p(:,index);
             q = q(:,index);
             w = w(:,index);
             totalw = sum(w,2);
             nrpoints = size(p,2);
         end
         %initialize(obj,p);
         obj.c = p;
         p = p-repmat(obj.c,1,nrpoints);% center p around its gravity
         q = q-repmat(obj.c,1,nrpoints);% center p around its gravity
         % align centers
         centerq = sum(w.*q,2)./totalw;
         centerp = sum(w.*p,2)./totalw;
         q = q-repmat(centerq,1,nrpoints);
         p = p-repmat(centerp,1,nrpoints);
%          % Initialization
%          mintheta=0.0001;
%          theta=pi/4;
%          transform=eye(3);
%          % iterate rotation
%          totalerror = calcep(q,p,w);
%          while theta>mintheta
%                M=calcmoment(q,p,w);
%                if (M==0), break; end
%                M=M/norm(M);
%                rotation = rodrigues(M*theta);
%                newp=rotation*p;
%                error=calcep(q,newp,w);
%                if (error<totalerror)
%                    p=newp;
%                    transform=rotation*transform;
%                    totalerror=error;
%                else
%                    theta=theta/2;
%                end 
%          end
         % rotation
         % Rotate
         [U,S,V] = svd(w.*q*p');
         transform = (V*sign(S)*U')';
         % putting it all together
         Ta=eye(4);
         Ta(1:3,4)=centerq;
         Tb=eye(4);
         Tb(1:3,4)=-centerp;
         R=eye(4);
         R(1:3,1:3)=transform;
         Tout=Ta*R*Tb;
         obj.Rotation = Tout(1:3,1:3);
         obj.Translation = Tout(1:3,4);
end

% function ep=calcep(q,p,w)
%          ep = sum(sum(w.*((p-q).^2)));
% end
% 
% function M=calcmoment(q,p,w)
%     M=sum(w.*cross(p,(q-p)),2);
% end


    

