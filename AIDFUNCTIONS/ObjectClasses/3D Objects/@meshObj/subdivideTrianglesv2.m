function varargout = subdivideTrianglesv2(obj,mode,val,varargin)
    if nargout == 1
           obj = clone(obj);
           obj.Visible = false;
           varargout{1} = obj;
    end
    [Vindex,Findex] = getVindexFindex(obj,varargin{:});
    if ~isempty(find(strcmp(varargin,'CheckFaces'), 1))
        if validChild(obj,obj.Map)
            in = inClusters(obj.Map,obj.Map.UV,0.01);% test 1
            vindex = find(in==0);
            findex = Vindex2Findex(obj,vindex);
            Findex = setdiff(Findex,findex);clear vindex findex;
            out = sameClusterFaces(obj.Map,obj.Faces(:,Findex));% test 2
            findex = find(out);
            Findex = intersect(Findex,findex);
%             good = goodUVFaces(obj.Map,obj.Faces(:,Findex),0.8);% test 3
%             findex = find(good);
%             Findex = intersect(Findex,findex);          
        end
    end
% switching mode    
    switch mode
        case 'runs'
            if (val<1), return; end
            for r=1:1:val
                Findex = performSubdivision(obj,Findex);
            end
        case 'size'
            if val<=0, return; end
            for r=1:1:50% safety maximum number of iterations
                face = obj.Faces(:,Findex);
                nrF = length(Findex);
                LOC =  zeros(3,nrF,3);AB = zeros(3,nrF);AC = zeros(3,nrF);
                for i=1:1:3
                    LOC(:,:,i) = reshape(obj.Vertices(i,face(:)),3,size(face,2));
                    AB(i,:) = LOC(1,:,i)-LOC(2,:,i);
                    AC(i,:) = LOC(1,:,i)-LOC(3,:,i);
                end
                areas = 0.5*sqrt(dot(AB,AB).*dot(AC,AC)-dot(AB,AC).^2);
                Aindex = find(areas>val);
                if isempty(Aindex), break; end% no triangles to subdivide
                Findex = Findex(Aindex);
                Findex = performSubdivision(obj,Findex);              
            end
            if(r==50), display('maximum number of iterations');end
        otherwise
            return;
    end
    updateChildren(obj,'Subdivide Triangles');
end

function Findex = performSubdivision(obj,Findex)
    %Findex = (1:1:size(obj.Faces,2));
    if isempty(Findex), return; end
    % one run of subdivision
    face = obj.Faces(:,Findex);
    
    
    
    
    
    face = obj.Faces(:,Findex);
    keepFindex = setdiff((1:1:size(obj.Faces,2)),Findex);
    vertex = obj.Vertices;
    n = size(vertex,2);
    % TRIANGLES
    i = [face(1,:) face(2,:) face(3,:) face(2,:) face(3,:) face(1,:)];
    j = [face(2,:) face(3,:) face(1,:) face(1,:) face(2,:) face(3,:)];
    I = find(i<j);
    i = i(I); j = j(I);
    [tmp,I] = unique(i + 1234567*j);
    i = i(I); j = j(I);
    ne = length(i); % number of edges
    s = n+(1:ne);
    
    A = sparse([i;j],[j;i],[s;s],n,n);
    try
        v12 = full( A( face(1,:) + (face(2,:)-1)*n ) );
        v23 = full( A( face(2,:) + (face(3,:)-1)*n ) );
        v31 = full( A( face(3,:) + (face(1,:)-1)*n ) );
    catch %#ok<CTCH>
        disp('Maximum number of Triangles obtained');
        Findex = [];
        return;
    end
        
    face = [   cat(1,face(1,:),v12,v31),...
                cat(1,face(2,:),v23,v12),...
                cat(1,face(3,:),v31,v23),...
                cat(1,v12,v23,v31)   ];
    
    % VERTICES
    obj.Vertices = [vertex, (vertex(:,i)+vertex(:,j))/2 ];
    clear vertex;
    obj.Faces = [obj.Faces(:,keepFindex) face];
    Findex = (1:size(face,2))+length(keepFindex);
    clear face;
    % GRADIENT
%     if ~isempty(obj.Gradient)
%        g = obj.Gradient;
%        g = [g, (g(:,i)+g(:,j))/2 ];
%        obj.Gradient = g;
%        clear g;
%     end
    % COLOR INFORMATION    
    if ~isempty(obj.TextureColor)
       rgb = obj.TextureColor;
       rgb = [rgb, (rgb(:,i)+rgb(:,j))/2 ];
       obj.TextureColor = rgb;
       clear rgb;
    end
    if ~isempty(obj.UV)
        uv = obj.UV;
        uv = [uv (uv(:,i)+uv(:,j))/2 ];
        obj.UV = uv;
        clear uv;
%         oldTex = obj.TextureColor;
%         getMapValues(obj,'Texture');
%         obj.TextureColor = replaceNanValues(obj,obj.TextureColor,'given',oldTex);
    end
end