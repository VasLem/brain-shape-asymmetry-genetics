function out = warp(obj,startUV,endUV,startC,endC)


    if nargout == 1
       obj1 = clone(obj1); 
       varargout{1} = obj1; 
    end
    nrC = length(obj1.UVClusters);
    if nrC == 1;% Simple warping of images
        cluster.Warped = warpImages(obj1,obj1.UV,obj2.UV);
    elseif nrC == 2;% complex warping of individual clusters
        cluster = warpClusters(obj1,obj2.UV);
        cluster = clusterInfluence(cluster);
    else
        return; % for now (change in future)
    end    
    % generating images from warped clusters
    Location = [];
    Texture = [];
    Indexed = [];
    Vertex = [];
    Gradient = [];
    for c=1:1:nrC
        Location = [Location cluster(c).Warped.UV]; %#ok<AGROW>
        Texture = [Texture cluster(c).Warped.Texture]; %#ok<AGROW>
        Indexed = [Indexed cluster(c).Warped.Indexed]; %#ok<AGROW>
        Vertex = [Vertex cluster(c).Warped.Vertex]; %#ok<AGROW>
        Gradient = [Gradient cluster(c).Warped.Gradient]; %#ok<AGROW>
        
    end   
    if~isempty(Texture),obj1.Texture = createImage(obj2.Frame,Location,Texture);clear Texture;end
    if~isempty(Indexed),obj1.Indexed = createImage(obj2.Frame,Location,Indexed);clear Texture;end
    if~isempty(Vertex),obj1.Vertex = createImage(obj2.Frame,Location,Vertex);clear Texture;end
    if~isempty(Gradient),obj1.Gradient = createImage(obj2.Frame,Location,Gradient);clear Texture;end
    obj1.UV = obj2.UV;
    obj1.Frame = obj2.Frame;
    obj1.UVClusters = obj2.UVClusters;
    obj1.Type = obj2.Type;
end

function cluster = clusterInfluence(cluster)
        nrC = length(cluster);
        FullIndex = (1:nrC);
        for k=1:1:nrC
            OtherIndex = setdiff(FullIndex,k);
            %figure;fastrbf_view(warped_cluster(k).uv,warped_cluster(k).Color.Value);
            for j=1:1:length(OtherIndex)
                eval = fastrbf_pointeval(cluster(OtherIndex(j)).ProbRBF,cluster(k).Warped.UV,'messages',0);
                %eval = fastrbf_pointeval(warped_cluster(other_index(j)).border.rbf,warped_cluster(k).uv,'messages',0);
                KeepIndex = find((cluster(k).Prob-eval.Value)>0);
                %keep_index = find((warped_cluster(k).border.Value-eval.Value)>0);
                clear eval;
                cluster(k).Warped.UV = cluster(k).Warped.UV(:,KeepIndex);
                cluster(k).Prob = cluster(k).Prob(KeepIndex);          
                if ~isempty(cluster(k).Warped.Texture)
                    cluster(k).Warped.Texture = cluster(k).Warped.Texture(:,KeepIndex);
                end
                if ~isempty(cluster(k).Warped.Indexed)
                    cluster(k).Warped.Indexed = cluster(k).Warped.Indexed(:,KeepIndex);
                end
                if ~isempty(cluster(k).Warped.Vertex)
                    cluster(k).Warped.Vertex = cluster(k).Warped.Vertex(:,KeepIndex);
                end
                if ~isempty(cluster(k).Warped.Gradient)
                    cluster(k).Warped.Gradient = cluster(k).Warped.Gradient(:,KeepIndex);
                end
                clear KeepIndex;
            end
        end
end

function map = createImage(Frame,Location,Value)
         nrV = size(Value,1);
         map = zeros(Frame.Dim(1),Frame.Dim(2),nrV);      
         for i=1:1:nrV
            map(:,:,i) = griddata(Location(1,:),Location(2,:),Value(i,:),Frame.X,Frame.Y,'linear');
         end
end