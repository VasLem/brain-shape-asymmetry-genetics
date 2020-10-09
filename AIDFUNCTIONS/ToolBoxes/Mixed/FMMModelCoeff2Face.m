function face = FMMModelCoeff2Face(model,coeff)
         vec = model.gem + model.eigvec*coeff'; 
         struc = vec2struc(vec,[],model.nr,model.Tri);
         face = meshObj(struc.mesh);
         try 
            face.TextureColor = struc.mesh.Color.Value; 
         catch
         end
         try 
            face.UserData.Properties = struc.properties; 
         catch
         end

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