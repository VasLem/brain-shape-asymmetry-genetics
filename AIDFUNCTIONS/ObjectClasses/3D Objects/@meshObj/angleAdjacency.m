function AA = angleAdjacency(obj)
% generates Adjancency matrix with angles between connected points
         face = obj.Faces;
         vertex = obj.Vertices;
         n = max(face(:));
         AA = sparse(n,n);
         vfring = computeVertexFaceRing(obj);
         for i=1:1:n
%              if mod(i,100)==0
%                  disp(num2str(i));
%              end
             for b = vfring{i}
                 bf = face(:,b);
                 % compute complementary vertices
                 if bf(1)==i
                    v = bf(2:3);
                 elseif bf(2)==i
                    v = bf([1 3]);
                 elseif bf(3)==i
                    v = bf(1:2);
                 else
                    error('Problem in face ring.');
                 end
                 j = v(1); k = v(2);
                 vi = vertex(:,i);
                 vj = vertex(:,j);
                 vk = vertex(:,k);
                 % angles
                 alpha = myangle(vk-vi,vk-vj);
                 beta = myangle(vj-vi,vj-vk);
                 % add weight
                 AA(i,j) = AA(i,j) + cot( alpha );
                 AA(i,k) = AA(i,k) + cot( beta );                 
             end
         end
         AA = (AA+AA')/2;
end

function beta = myangle(u,v)
    du = sqrt( sum(u.^2) );
    dv = sqrt( sum(v.^2) );
    du = max(du,eps); dv = max(dv,eps);
    beta = acos( sum(u.*v) / (du*dv) );
end