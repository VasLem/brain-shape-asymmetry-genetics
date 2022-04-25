function importWrl(obj,filename,varargin)
    % changing directory
         Input = find(strcmp(varargin,'Path'));
         if ~isempty(Input)
            cd(varargin{Input+1});
         end
     % Getting type
         Input = find(strcmp(varargin,'Type'));
         if isempty(Input)
            type = 'Wrl';
         else
            type = varargin{Input+1};
         end
         [vertex,face,uv,uvface] = read_wrl(filename);
         
      % Storing information in object
          %obj.ColorMode = 'Single';
          obj.Vertices = vertex;
          obj.Faces = face;
          %obj.ColorMode = 'Single';
      % storing kind of wavefront in Userdatas
          obj.UserData = lower(type);
      % reading texture file
          if isempty(uv), return; end % there is no texture image to look for
          name = filename(1:end-4);
          tfile = dir([name '.bmp']);
          if isempty(tfile), return; end% texture image does not exists
          im = double(imread(tfile.name))./255;
          factor = round(numel(im)/4000000);
          if factor > 1, im = im(1:factor:end,1:factor:end,:); end
          uv(2,:) = 1-uv(2,:);
          obj.UV = uv;
          obj.TextureMap = im;
          obj.ColorMode = 'Texture';
end


function [vertex, face, uv, uvface] = read_wrl(filename)

% read_wrl - load a mesh from a VRML file
%
%   [vertex, face] = read_wrl(filename);
%
%   Copyright (c) 2004 Gabriel Peyré
    uv = [];
    uvface = [];
    fid = fopen(filename,'r');
    data = textscan(fid,'%[^\n\r]');data = data{1};
    fclose(fid);
    CloseBrackets = strmatch(']',data);
    OpenPoints = strmatch('point [',data);
    startvertex = OpenPoints(1)+1;
    tmp = find(CloseBrackets>startvertex);
    endvertex = CloseBrackets(tmp(1))-1;
    % reading vertices
    index = (startvertex:endvertex);
    vertex = nan*zeros(3,length(index));
    for i=1:1:length(index)
        C = textscan(data{index(i)},'%f');
        vertex(:,i) = C{1};
    end
    clear index;
    startface = strmatch('coordIndex [',data)+1;
    tmp = find(CloseBrackets>startface);
    endface = CloseBrackets(tmp(1))-1;
    index = (startface:endface);
    TRI = nan*zeros(3,length(index));
    QUAD = nan*zeros(4,length(index));
    for i=1:1:length(index)
        str = data{index(i)};
        str = cleanStringFrom(str,',');
        C = textscan(str,'%f');
        switch length(C{1})
            case 4
                TRI(:,i) = C{1}(1:end-1);
            case 5
                QUAD(:,i) = C{1}(1:end-1);
        end
    end
    clear index;
    TRI(:,find(isnan(TRI(1,:)))) = [];
    QUAD(:,find(isnan(QUAD(1,:)))) = [];
    face = TRI;
    if ~isempty(QUAD)
        Tri1 = zeros(3,size(QUAD,2));
        Tri2 = zeros(3,size(QUAD,2));
        for k=1:1:size(QUAD,2)
            Tri1(:,k) = QUAD(1:3,k);
            Tri2(1,k) = QUAD(1,k);
            Tri2(2:end,k) = QUAD(3:end,k);
        end
       face = [face Tri1 Tri2];
    end
    face = face+1;
    clear TRI QUAD;
    if length(OpenPoints)>1
       % then we have texture coordinates!!!
       startuv = OpenPoints(2)+1;
       tmp = find(CloseBrackets>startuv);
       enduv = CloseBrackets(tmp(1))-1;
       index = (startuv:enduv);
       uv = nan*zeros(2,length(index));
       for i=1:1:length(index)
           C = textscan(data{index(i)},'%f');
           uv(:,i) = C{1};
       end
       clear index;
       startuvface = strmatch('texCoordIndex [',data)+1;
       tmp = find(CloseBrackets>startuvface);
       enduvface = CloseBrackets(tmp(1))-1;
       index = (startuvface:enduvface);
       TRI = nan*zeros(3,length(index));
       QUAD = nan*zeros(4,length(index));
       for i=1:1:length(index)
           str = data{index(i)};
           str = cleanStringFrom(str,',');
           C = textscan(str,'%f');
           switch length(C{1})
              case 4
                   TRI(:,i) = C{1}(1:end-1);
              case 5
                   QUAD(:,i) = C{1}(1:end-1);
           end
       end
       clear index;
       TRI(:,find(isnan(TRI(1,:)))) = [];
       QUAD(:,find(isnan(QUAD(1,:)))) = [];
       uvface = TRI;
       if ~isempty(QUAD)
          Tri1 = zeros(3,size(QUAD,2));
          Tri2 = zeros(3,size(QUAD,2));
          for k=1:1:size(QUAD,2)
              Tri1(:,k) = QUAD(1:3,k);
              Tri2(1,k) = QUAD(1,k);
              Tri2(2:end,k) = QUAD(3:end,k);
          end
          uvface = [uvface Tri1 Tri2];
       end
       uvface = uvface+1;
       
       [uniqF, uniqI] = unique(face(:));
       tmpuvface = uvface(:);
       tmpuv = zeros(2,size(vertex,2));
       tmpuv(:,uniqF) = uv(:,tmpuvface(uniqI));
       uv = tmpuv;
    end
end


%     fid = fopen(filename,'r');
%     data = textscan(fid,'%[^\n\r]');data = data{1};
%     fclose(fid);
%     CloseBrackets = strmatch(']',data);
%     startvertex = strmatch('point [',data)+1;
%     % reading vertices
%     index = (startvertex:CloseBrackets(1)-1);
%     vertex = nan*zeros(3,length(index));
%     for i=1:1:length(index)
%         C = textscan(data{index(i)},'%f');
%         vertex(:,i) = C{1};
%     end
%     clear index;
%     startface = strmatch('coordIndex [',data)+1;
%     index = (startface:CloseBrackets(2)-1);
%     TRI = nan*zeros(3,length(index));
%     QUAD = nan*zeros(4,length(index));
%     for i=1:1:length(index)
%         str = data{index(i)};
%         str = cleanStringFrom(str,',');
%         C = textscan(str,'%f');
%         switch length(C{1})
%             case 4
%                 TRI(:,i) = C{1}(1:end-1);
%             case 5
%                 QUAD(:,i) = C{1}(1:end-1);
%         end
%     end
%     TRI(:,find(isnan(TRI(1,:)))) = [];
%     QUAD(:,find(isnan(QUAD(1,:)))) = [];
%     face = TRI;
%     if ~isempty(QUAD)
%         Tri1 = zeros(3,size(QUAD,2));
%         Tri2 = zeros(3,size(QUAD,2));
%         for k=1:1:size(QUAD,2)
%             Tri1(:,k) = QUAD(1:3,k);
%             Tri2(1,k) = QUAD(1,k);
%             Tri2(2:end,k) = QUAD(3:end,k);
%         end
%        face = [face Tri1 Tri2];
%     end
%     face = face+1;
