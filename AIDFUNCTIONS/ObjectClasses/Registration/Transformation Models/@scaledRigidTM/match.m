function out = match(obj,q,p,w)
% find the rigid transform from p to q
         if nargout == 1
            obj = clone(obj);
            out = obj;
         end
         p = TM.getPoints(p);
         q = TM.getPoints(q);
         nbpts = size(p,2);
         if nargin<4
            w  = ones(1,nbpts);
         end 
         totalw = sum(w);
         index = find(w);% find points not having weights == 0;
         if length(index)<nbpts
             p = p(:,index);
             q = q(:,index);
             w = w(:,index);
             totalw = sum(w,2);
             nbpts = size(p,2);
         end
         rw = repmat(w,3,1);
         % align centers
         centerq = sum(rw.*q,2)./repmat(totalw,3,1);
         centerp = sum(rw.*p,2)./repmat(totalw,3,1);
         q = q-repmat(centerq,1,nbpts);
         p = p-repmat(centerp,1,nbpts);
         % Scale both by making mean length 1
         lengths = sqrt(sum((rw.*p).^2,1));
         meanLength1 = sum(lengths,2)/totalw;
         p = p./meanLength1;
         lengths = sqrt(sum((rw.*q).^2,1));
         meanLength2 = sum(lengths,2)/totalw;
         q = q./meanLength2;
         scalefactor = meanLength2/meanLength1;
         % Rotate
         [U,S,V] = svd(rw.*q*p');
         H = V*sign(S)*U';
         H = H';
         % putting it all together
         obj.Scale=scalefactor;
         obj.Rotation = H;
         transform = obj.Scale*H;
         Ta=eye(4);
         Ta(1:3,4)=centerq;%+obj.c;
         Tb=eye(4);
         Tb(1:3,4)=-centerp;%-obj.c;
         R=eye(4);
         R(1:3,1:3)=transform;
         Tout=Ta*R*Tb;
         obj.Translation = Tout(1:3,4)/obj.Scale;
end
    



% function out = match(obj,q,p,w)
% % find the rigid transform from p to q
%          if nargout == 1
%             obj = clone(obj);
%             out = obj;
%          end
%          p = TM.getPoints(p);
%          q = TM.getPoints(q);
%          nrpoints = size(p,2);
%          if nargin<4
%             w  = ones(3,nrpoints);
%          else
%             if size(w,1)==1, w = repmat(w,3,1);end 
%          end
%          totalw = sum(w,2);
%          index = find(sum(w));% find points not having weights == 0;
%          if length(index)<nrpoints
%              p = p(:,index);
%              q = q(:,index);
%              w = w(:,index);
%              totalw = sum(w,2);
%              nrpoints = size(p,2);
%          end
%          %initialize(obj,p);
% %          obj.c = p;
% %          p = p-repmat(obj.c,1,nrpoints);% center p around its gravity
% %          q = q-repmat(obj.c,1,nrpoints);% center p around its gravity
%          % align centers
%          centerq = sum(w.*q,2)./totalw;
%          centerp = sum(w.*p,2)./totalw;
%          q = q-repmat(centerq,1,nrpoints);
%          p = p-repmat(centerp,1,nrpoints);
%          % Initialization
%          %iterate rotation
%           iteration=0;
%           endit=0;
%           transform=eye(3);
%           while (endit<3)
%             iteration=iteration+1;
%             Za=p(1,:)*(q(1,:))'+p(2,:)*(q(2,:))';
%             Zb=p(2,:)*(q(1,:))'-p(1,:)*(q(2,:))';
%             normZ=sqrt(Za*Za+Zb*Zb);
%             Za=Za/normZ;
%             Zb=Zb/normZ;
%             if (Za<1)
%                 endit=0;
%                 rotation=eye(3);
%                 rotation(1,1)=Za;
%                 rotation(1,2)=Zb;
%                 rotation(2,1)=-Zb;
%                 rotation(2,2)=Za;
%                 p=rotation*p;
%                 transform=rotation*transform;
%             else
%                 endit=endit+1;
%             end
%             %2
%             Ya=p(2,:)*(q(2,:))'+p(3,:)*(q(3,:))';
%             Yb=p(3,:)*(q(2,:))'-p(2,:)*(q(3,:))';
%             normY=sqrt(Ya*Ya+Yb*Yb);
%             Ya=Ya/normY;
%             Yb=Yb/normY;
%             if (Ya<1)
%                 endit=0;
%                 rotation=eye(3);
%                 rotation(2,2)=Ya;
%                 rotation(2,3)=Yb;
%                 rotation(3,2)=-Yb;
%                 rotation(3,3)=Ya;
%                 p=rotation*p;
%                 transform=rotation*transform;
%             else
%                 endit=endit+1;
%             end
%             %3
%             Xa=p(3,:)*(q(3,:))'+p(1,:)*(q(1,:))';
%             Xb=p(1,:)*(q(3,:))'-p(3,:)*(q(1,:))';
%             normX=sqrt(Xa*Xa+Xb*Xb);
%             Xa=Xa/normX;
%             Xb=Xb/normX;
%             if (Xa<1)
%                 endit=0;
%                 rotation=eye(3);
%                 rotation(3,3)=Xa;
%                 rotation(3,1)=Xb;
%                 rotation(1,3)=-Xb;
%                 rotation(1,1)=Xa;
%                 p=rotation*p;
%                 transform=rotation*transform;
%             else
%                 endit=endit+1;
%             end
%           end
%         %determine scale
%         In=sum(sum(p.*p));
%         obj.Scale=(p(1,:)*(q(1,:))'+p(2,:)*(q(2,:))'+p(3,:)*(q(3,:))')/In;
%         obj.Rotation = transform;
%         transform = obj.Scale*transform;
%         Ta=eye(4);
%         Ta(1:3,4)=centerq;%+obj.c;
%         Tb=eye(4);
%         Tb(1:3,4)=-centerp;%-obj.c;
%         R=eye(4);
%         R(1:3,1:3)=transform;
%         %centerp
%         Tout=Ta*R*Tb;
%         %obj.Rotation = Tout(1:3,1:3);
%         obj.Translation = Tout(1:3,4)/obj.Scale;
%         %obj.Translation = Tout(1:3,4);
% end
%     
% 
