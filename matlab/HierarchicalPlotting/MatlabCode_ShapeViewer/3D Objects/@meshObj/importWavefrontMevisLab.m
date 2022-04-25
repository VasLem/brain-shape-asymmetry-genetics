function importWavefrontMevisLab(obj,filename,varargin)
    % changing directory
         Input = find(strcmp(varargin,'Path'));
         if ~isempty(Input)
            cd(varargin{Input+1});
         end
     % Getting type
         Input = find(strcmp(varargin,'Type'));
         if isempty(Input)
            type = 'Wavefront';
         else
            type = varargin{Input+1};
         end
     % quick inreading of vertices, normals and texturecoordinates;
         fid = fopen(filename,'r');
         names = textscan(fid,'%s%*[^\n]');names = names{1};
         fclose(fid);
         fid = fopen(filename,'r');
         data = textscan(fid,'%[^\n\r]');
         data = data{1};
         fclose(fid);
     % reading Vertices
         %disp('reading vertices');
         list = strcmp('v',names);index = find(list==1);
         if isempty(index), return; end %nothing in the file
         % small test to check whether rgb values are given
         Lyn = data{index(1)}; tmp = sscanf(Lyn(2:end),'%f');
         Location = nan*zeros(size(tmp,1),length(index));
         for i=1:1:length(index)
             Lyn = data{index(i)};
             Location(:,i) = sscanf(Lyn(2:end),'%f');
         end
         vertex = Location(1:3,:);
         if size(Location,1)>3 
             texturecolor = Location(4:end,:)./255; 
         else
             texturecolor = [];
         end
         clear list index Location;
     % reading normals 
         % not required because computed in meshOBj class on the fly    
     % reading texture coordinates
         %disp('reading texture coordinates');
         list = strcmp('vt',names);index = find(list==1);
         if ~isempty(index)
            uv = nan*zeros(2,length(index));
            for i=1:1:length(index)
                Lyn = data{index(i)};
                uv(:,i) = sscanf(Lyn(3:end),'%f');
            end
         else
             uv = [];
         end  
         clear list index;
     % reading faces
        %disp('reading faces');
        list = strcmp('f',names);index = find(list==1);
        if ~isempty(index)
            F3 = nan*zeros(3,length(index));F3t = F3;
            F4 = nan*zeros(4,length(index));F4t = F4;
            for i=1:1:length(index)
                Lyn = data{index(i)};
                Lyn=deblank(Lyn(3:end));
                nvrts=length(findstr(Lyn,' ')); % nr of vertices
                nslash=length(findstr(Lyn,'/')); % nr of slashes
                %keyboard;
                switch nvrts
                    case 3 % Triangles
                        switch nslash
                            case 0 % vertex
                                F3(:,i) = sscanf(Lyn,'%f');
                            case 3 % vertex\texture
                                f=sscanf(Lyn,'%f/%f');
                                F3t(:,i) = f([2 4 6]);
                                F3(:,i) = f([1 3 5]);
                            case 6 % vertex\texture\normal
                                f=sscanf(Lyn,'%f/%f/%f');
                                %keyboard;
                                try
                                 F3t(:,i) = f([2 5 8]);
                                 F3(:,i) = f([1 4 7]);
                                catch
                                 f=sscanf(Lyn,'%f//%f');
                                 %keyboard;
                                 %F3t(:,i) = f([2 5 8]);
                                 F3(:,i) = f([1 3 6]);
                                end
                            otherwise
                        end
                    case 4 % Quadruples
                        switch nslash
                            case 0 % vertex
%                                 F4(:,i) = sscanf(Lyn,'%f');
                            case 4 % vertex\texture
                                f=sscanf(Lyn,'%f/%f');
                                F4t(:,i) = f([2 4 6 8]);
                                F4(:,i) = f([1 3 5 7]);
                            case 8 % vertex\texture\normal
                                f=sscanf(Lyn,'%f/%f/%f');
                                F4t(:,i) = f([2 5 8 11]);
                                F4(:,i) = f([1 4 7 10]);
                            otherwise
                        end
                    otherwise
                end
            end
            % trimming the face lists, dropping the nan values
            index = find(~isnan(F3(1,:)));F3 = F3(:,index);
            index = find(~isnan(F3t(1,:)));F3t = F3t(:,index);
            index = find(~isnan(F4(1,:)));F4 = F4(:,index);
            index = find(~isnan(F4t(1,:)));F4t = F4t(:,index);
            %converting quads to triangles
            if ~isempty(F4)
                Tri1 = zeros(3,size(F4,2));
                Tri2 = zeros(3,size(F4,2));
                        for k=1:1:size(F4,2)
                            Tri1(:,k) = F4(1:3,k);
                            Tri2(1,k) = F4(1,k);
                            Tri2(2:end,k) = F4(3:end,k);
                        end
                F3 = [F3,Tri1,Tri2];
            end
            if ~isempty(F4t)
                Tri1 = zeros(3,size(F4t,2));
                Tri2 = zeros(3,size(F4t,2));
                        for k=1:1:size(F4t,2)
                            Tri1(:,k) = F4t(1:3,k);
                            Tri2(1,k) = F4t(1,k);
                            Tri2(2:end,k) = F4t(3:end,k);
                        end
                F3t = [F3t,Tri1,Tri2];
            end
        else
            F3 = [];F3t = F3;
            F4 = [];F4t = F4;
        end
     % Deleting non mesh points
        in_mesh = unique(F3(:))';
        not_nan = find(~isnan(vertex(1,:)));
        good = intersect(in_mesh,not_nan);
        %keyboard;
     % adjusting triangles
        [tmp,LOC] = ismember(F3,good);
        [i,j] = find(tmp==0);
        tmpindex = (1:size(F3,2));
        tmpindex = setdiff(tmpindex,j);
        F3 = LOC(:,tmpindex);
     % adjusting texture coord information
        if ~isempty(F3t), F3t = F3t(:,tmpindex); end
     % adjusting vertex information    
        vertex = vertex(:,good);
     %adjusting rgb information
        if ~isempty(texturecolor), texturecolor = texturecolor(:,good); end
     % rearrangig texture coordinates
        if ~isempty(uv)           
           if ~isempty(F3t)
              [tmp,I,J] = unique(F3(:));
              T_index = F3t(I);
              uv = uv(:,T_index);    
           else
              uv = uv(:,good);
           end
        end
      % Storing information in object
          obj.ColorMode = 'Single';
          obj.Vertices = vertex;
          obj.Faces = F3;
          if ~isempty(texturecolor)
              obj.TextureColor = texturecolor;
              obj.ColorMode = 'Texture';
          end
      % storing kind of wavefront in Userdatas
          obj.UserData = lower(type);
      % reading texture file
          if isempty(uv), return; end % there is no texture image to look for
          name = filename(1:end-4);
          if ~strcmp(type,'Eyetronics')
            tfile = dir([name '.bmp']);
          else
            tfile = dir([name '_tex.bmp']);  
          end
%           if isempty(tfile), 
%              tfile = dir([name '.gif']);
%           end
          if isempty(tfile),return; end% texture image does not exists
          im = double(imread(tfile.name))./255;
          factor = round(numel(im)/4000000);
          if factor > 1, im = im(1:factor:end,1:factor:end,:); end
          uv(2,:) = 1-uv(2,:);
          %uv(1,:) = 1-uv(1,:);
          obj.UV = uv;
          obj.TextureMap = im;
          obj.ColorMode = 'Texture';
end