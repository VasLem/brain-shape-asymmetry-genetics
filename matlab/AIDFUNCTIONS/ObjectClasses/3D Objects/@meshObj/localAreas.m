function out = localAreas(obj)
          face = obj.Faces;
          nrF = size(face,2);
          LOC =  zeros(3,nrF,3);AB = zeros(3,nrF);AC = zeros(3,nrF);
          for i=1:1:3
                   LOC(:,:,i) = reshape(obj.Vertices(i,face(:)),3,size(face,2));
                   AB(i,:) = LOC(1,:,i)-LOC(2,:,i);
                   AC(i,:) = LOC(1,:,i)-LOC(3,:,i);
          end
          areas = 0.5*sqrt(dot(AB,AB).*dot(AC,AC)-dot(AB,AC).^2);
          out = zeros(1,obj.nrV);
          parfor i=1:obj.nrV
                 tmp = ismember(obj.Faces,i);
                 [~,j] = find(tmp==1);
                 out(i) = mean(areas(j));
          end
end