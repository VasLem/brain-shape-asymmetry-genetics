function [Data,RefScan] = FMMModel2Data(model)
         if isfield(model,'name'),Data.Names = model.name;end
         if isfield(model,'property'), Data.Properties.Names = model.property.name;end
         nrsamples = size(model.bcoeff,1);
         if model.nr(1)>0, Data.Shape = zeros(model.nr(1),nrsamples); end
         if model.nr(2)>0, Data.TextureColor = uint8(zeros(model.nr(2),nrsamples)); end
         if model.nr(3)>0, Data.Properties.Values = zeros(model.nr(3),nrsamples); end
         for i=1:1:size(model.bcoeff,1)
             vec = model.gem + model.eigvec*model.bcoeff(i,:)';
             struc = vec2struc(vec,[],model.nr,model.Tri);
             if model.nr(1)>0,Data.Shape(:,i) = struc.mesh.Location(:);end
             if model.nr(2)>0,Data.TextureColor(:,i) = uint8(round(struc.mesh.Color.Value(:)*255));end
             if model.nr(2)>0,Data.Properties.Value(:,i) = struc.properties;end
             clear vec struc;
         end
         RefScan = meshObj(model.gem_face.mesh);
         RefScan.TextureColor = model.gem_face.mesh.Color.Value;
end
function [structure] = vec2struc(vector,bcoeff,nr,Tri)
    if size(vector,2) > size(vector,1)
        vector = vector';
    end
    if size(bcoeff,2) > size(bcoeff,1)
        bcoeff = bcoeff';
    end

    structure.vector = vector;
    if ~isempty(bcoeff)
        structure.bcoeff = bcoeff;
    end   
    if ~(nr(3)==0)
       % Then we have properties
       structure.properties = vector(end-nr(3)+1:end);
       vector = vector(1:end-nr(3));
    end
    if ~(nr(2)==0)
        % we have texture information
        color_vector = vector(end-nr(2)+1:end);
        vector = vector(1:end-nr(2));
        color_vector = reshape(color_vector,3,size(color_vector,1)/3);
        structure.mesh.Color.Value = [color_vector(1,:);color_vector(2,:);color_vector(3,:)];
    end
    if ~(nr(1)==0)
        % we have shape information
        mesh_vector = vector(1:nr(1));
        structure.mesh.Location = reshape(mesh_vector,3,size(mesh_vector,1)/3);
        structure.mesh.Tri = Tri;
    end
end