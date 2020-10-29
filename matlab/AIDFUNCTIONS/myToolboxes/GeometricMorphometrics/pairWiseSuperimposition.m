function Floating = pairWiseSuperimposition(Target,Floating,scale,kappa,index)
         if nargin<4, kappa = +inf;end% kappa = +inf -> simple Least Squares Superimposition
         if nargin<5, index = 1:Target.nVertices;end
         if nargin<3, scale = false; end
         if nargout==1, Floating = clone(Floating); end% clone outcome, not to change the input
         if ~(Target.nVertices==Floating.nVertices), error('A superimposition cannot be done between shapes that do not share the same amount of vertices'); end       
         display = false;   
         if kappa==+inf % simple unweighted procrustes superimposition
             % compute transformation
              [~,~,transform] = procrustes(Target.Vertices(index,:),Floating.Vertices(index,:),'Scaling',scale,'Reflection',false);
             % evaluate transformation 
               Floating.Vertices = transform.b*Floating.Vertices*transform.T + repmat(transform.c(1,:),Floating.nVertices,1);
            return;% done and return
         end
         % initialize inlier beliefs
         winit = zeros(1,Floating.nVertices);
         winit(index)= 1;% activate only points in the index
         w = winit;sumw = sum(w);
         if display
            Floating.VertexValue = w;
            v = viewer(Floating); viewer(Target,v);
            Target.ViewMode = 'points';
            Target.SingleColor = [0.8 0.8 0.8];
            Target.ColorMode = 'single';
            Floating.ViewMode = 'Wireframe';
            Floating.ColorMode = 'Indexed';
            pause(1);
         end
         for iter = 1:1:10% perform ten times
            % get transformation
             T = getTransformation(Target,Floating,scale,w);
            % evaluate transformation
             evalTransformation(T,Floating);
            % get remaining distances
             dist = sqrt(sum((Floating.Vertices-Target.Vertices).^2,2))';
            % update sigma of inlier gaussian pdf
             sigma = sqrt(sum(w.*dist.^2)/sumw);
            % update Level outlier uniform pdf  
             L = (1/sqrt(((2*pi)^2)*det(sigma)))*exp(-0.5*kappa^2);
            % update inlier beliefs
             ip = normpdf(dist,0,sigma);
            % update outlier beliefs 
             op = repmat(L,1,Floating.nVertices);
            % update weights 
             w = ip./(ip+op);
             w = w.*winit;sumw = sum(w);% incorporate deterministic index by dot multiplication with winit
             if display, Floating.VertexValue = w;pause(0.1); end
         end
         Floating.UserData.InlierBeliefs = w;
end

function T = getTransformation(Target,Floating,scale,w)
         q = Target.Vertices';
         p = Floating.Vertices';
         nbpts = size(p,2);
         index = find(w);% find points not having weights == 0;
         if length(index)<nbpts
             p = p(:,index);
             q = q(:,index);
             w = w(:,index);
         end
         totalw = sum(w,2);
         rw = repmat(w,3,1);
         nbpts = size(p,2);
         % align centers
         centerq = sum(rw.*q,2)./repmat(totalw,3,1);
         centerp = sum(rw.*p,2)./repmat(totalw,3,1);
         q = q-repmat(centerq,1,nbpts);
         p = p-repmat(centerp,1,nbpts);
         % Scale both by making mean length 1
         if scale
             lengths = sqrt(sum((rw.*p).^2,1));
             meanLength1 = sum(lengths,2)/totalw;
             p = p./meanLength1;
             lengths = sqrt(sum((rw.*q).^2,1));
             meanLength2 = sum(lengths,2)/totalw;
             q = q./meanLength2;
             scalefactor = meanLength2/meanLength1;
         else
             scalefactor = 1;% keep scale fixed
         end
         % Rotate
         [U,S,V] = svd(rw.*q*p');
         H = V*sign(S)*U';
         H = H';
         % putting it all together
         T.Scale=scalefactor;
         T.Rotation = H;
         transform = T.Scale*H;
         Ta=eye(4);
         Ta(1:3,4)=centerq;
         Tb=eye(4);
         Tb(1:3,4)=-centerp;
         R=eye(4);
         R(1:3,1:3)=transform;
         Tout=Ta*R*Tb;
         T.Translation = Tout(1:3,4)/T.Scale;
end

function evalTransformation(T,Floating)
         Floating.Vertices = (T.Scale*(T.Rotation*(Floating.Vertices')+repmat(T.Translation,1,Floating.nVertices)))';
end