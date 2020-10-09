function [out,totalerror]=rigidMatch(a,b,varargin)
% match(IN1,IN2) computes the optimal 4x4 homogenous transformation matrix M,
% so that M*IN2 approximates IN1 in a least squares sense
%
%                   _           _
%                  | R3x3   T3x1 |
%              M = |_  0      1 _|
%
%
% IN1 and IN2 are matrices containing a 3D point in each column, row1 contains the 
% x-, row2 the y- and row3 the z-coordinates of the points to be matched
%
% A 3D point in column x of IN1 will be matched with the 3D point in column x of IN2
%
% CODE: IN1(:,x) will be matched with IN2(:,x)
%
% The number of points that are matched is equal to the minimum value of the number 
% of points in IN1 and IN2
%
% CODE: nrpoints=min(size(IN1,2),size(IN2,2))
%
% IMPORTANT: For the computation of M*IN2, IN2 has to be converted to its homogeneous
%            equivalent by adding a row of ones. After the transformation has been executed,
%            this additional row has to be stripped from the resulting matrix
%            
%            CODE:        IN2HOM=[IN2;ones(1,nrpoints];
%                         NEWIN2HOM=M*IN2HOM;
%                         NEWIN2=IN2HOM(3,:);
%
% If the number of points to be matched equals zero, the resulting matrix M will be 
% the 4x4 Identity matrix

Input = find(strcmp(varargin, 'Weights'));
if ~isempty(Input)
    [out,totalerror]=rigidMatchW(a,b,varargin{Input+1});
    return;
end


nrpoints=min(size(a,2),size(b,2));

if (nrpoints>0)
    mintheta=0.0001;
    copya=a;
    copyb=b;

    theta=pi/4;
    transform=eye(3);

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
    totalerror=calcep(a,b);

    while (theta>mintheta)
        M=calcmoment(a,b);
        if (M==0)
            break;
        end
        M=M/norm(M);
        rotation=rodrigues(M*theta);
        c=rotation*b;
    
        error=calcep(a,c);
        if (error<totalerror)
            b=c;
            transform=rotation*transform;
            totalerror=error;
        else
            theta=theta/2;
        end
    end

    Ta=eye(4);
    Ta(1:3,4)=centera;
    Tb=eye(4);
    Tb(1:3,4)=-centerb;
    R=eye(4);
    R(1:3,1:3)=transform;

    out=Ta*R*Tb;
else
    out= eye(4);
end

% figure(1);
% hold on;
% plotpoints(copya,'rx');
% plotpoints(copyb,'gx');
% plotpoints(out*[copyb;ones(1,nrpoints)],'bo');
% hold off;
return;