function exportOff(obj,filename,varargin)
    % changing directory
         Input = find(strcmp(varargin,'Path'));
         if ~isempty(Input)
            cd(varargin{Input+1});
         end
     % Renormalize?
         if ~isempty(find(strcmp(varargin,'Renormalize'), 1))
            renormalize = 1;
         else
            renormalize = 0;
         end
         vertex = obj.Vertices;
         face = obj.Faces;
         if size(vertex,2)~=3
            vertex=vertex';
         end
         if size(vertex,2)~=3
            error('vertex does not have the correct format.');
         end

         if renormalize==1
             m = mean(vertex);
             s = std(vertex);
             for i=1:3
                 vertex(:,i) = (vertex(:,i)-m(i))/s(i);
             end
         end

         if size(face,2)~=3
             face=face';
         end
         if size(face,2)~=3
             error('face does not have the correct format.');
         end

         fid = fopen(filename,'wt');
         if( fid==-1 )
             error('Can''t open the file.');
             return;
         end

         % header
         fprintf(fid, 'OFF\n');
         fprintf(fid, '%d %d 0\n', size(vertex,1), size(face,1));

         % write the points & faces
         fprintf(fid, '%f %f %f\n', vertex');
         fprintf(fid, '3 %d %d %d\n', face'-1);

         fclose(fid);
     
end