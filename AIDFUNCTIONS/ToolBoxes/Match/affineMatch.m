function out=affineMatch(a,b)
% rotatie translatie en schaalfactor

copya=a;
copyb=b;

transform=eye(3);
nrpoints=size(a,2);
c=zeros(size(a));

%allign centers
centera=zeros(3,1);
centerb=zeros(3,1);
for i=1:nrpoints
    centera=centera+a(:,i);
    centerb=centerb+b(:,i);
end
centera=centera/nrpoints;
centerb=centerb/nrpoints;

for i=1:nrpoints
    a(:,i)=a(:,i)-centera;
    b(:,i)=b(:,i)-centerb;
end

%iterate rotation
iteration=0;

endit=0;

while (endit<3)
    %calcep(a,b)

    iteration=iteration+1;
    
    Za=b(1,:)*(a(1,:))'+b(2,:)*(a(2,:))';
    Zb=b(2,:)*(a(1,:))'-b(1,:)*(a(2,:))';
    normZ=sqrt(Za*Za+Zb*Zb);
    Za=Za/normZ;
    Zb=Zb/normZ;
    
    if (Za<1)
        endit=0;
    
        rotation=eye(3);
        rotation(1,1)=Za;
        rotation(1,2)=Zb;
        rotation(2,1)=-Zb;
        rotation(2,2)=Za;
        b=rotation*b;
   
        transform=rotation*transform;
    else
        endit=endit+1;
    end
    
    %2
    Ya=b(2,:)*(a(2,:))'+b(3,:)*(a(3,:))';
    Yb=b(3,:)*(a(2,:))'-b(2,:)*(a(3,:))';
    normY=sqrt(Ya*Ya+Yb*Yb);
    Ya=Ya/normY;
    Yb=Yb/normY;
    
    if (Ya<1)
        endit=0;
       
        rotation=eye(3);
        rotation(2,2)=Ya;
        rotation(2,3)=Yb;
        rotation(3,2)=-Yb;
        rotation(3,3)=Ya;
    
        b=rotation*b;
    
        transform=rotation*transform;
    else
        endit=endit+1;
    end
    
    %3
    Xa=b(3,:)*(a(3,:))'+b(1,:)*(a(1,:))';
    Xb=b(1,:)*(a(3,:))'-b(3,:)*(a(1,:))';
    normX=sqrt(Xa*Xa+Xb*Xb);
    Xa=Xa/normX;
    Xb=Xb/normX;
    
    if (Xa<1)
        endit=0;
        
        rotation=eye(3);
        rotation(3,3)=Xa;
        rotation(3,1)=Xb;
        rotation(1,3)=-Xb;
        rotation(1,1)=Xa;
    
        b=rotation*b;
    
        transform=rotation*transform;
    else
        endit=endit+1;
    end
end

%determine scale
In=sum(sum(b.*b));
scaletransform=(b(1,:)*(a(1,:))'+b(2,:)*(a(2,:))'+b(3,:)*(a(3,:))')/In;
transform=scaletransform*transform;


%concatenate result
Ta=eye(4);
Ta(1:3,4)=centera;
Tb=eye(4);
Tb(1:3,4)=-centerb;
R=eye(4);
R(1:3,1:3)=transform;

out=Ta*R*Tb;

% figure(1);
% hold on;
% plotpoints(copya,'rx');
% plotpoints(copyb,'gx');
% plotpoints(out*[copyb;ones(1,nrpoints)],'bo');
% 
% iteration
end