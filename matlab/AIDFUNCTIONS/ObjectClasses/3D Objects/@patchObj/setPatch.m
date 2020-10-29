function setPatch(obj,action)
if isempty(obj.ph)||~ishandle(obj.ph)||~obj.UpdatePatch, return; end
% if isempty(obj.ph)||~ishandle(obj.ph), return; end
switch action
    case 'All'
        setPatch(obj,'Vertices');
        setPatch(obj,'ColorMode');
        setPatch(obj,'ViewMode');
        setPatch(obj,'Faces');
%         setPatch(obj,'Gradient');
        setPatch(obj,'Viewer');
        setPatch(obj,'Transparancy');
        setPatch(obj,'Material');
        setPatch(obj,'LightMode');
        setPatch(obj,'MarkerSize');
        setPatch(obj,'Visible');
        return;
    case 'Shape'
        setPatch(obj,'Vertices');
        setPatch(obj,'Faces');
%         setPatch(obj,'Gradient');
    case 'VerticesColor'
        setPatch(obj,'Vertices');
        setPatch(obj,'Faces');
        setPatch(obj,'ColorMode');
        return;
    case 'Axes'
        prop = {'Parent'};val{1} = obj.Axes;
    case 'Vertices'
        %if isempty(obj.Vertices), return; end
        prop = {'Vertices'}; val{1} = obj.Vertices';
    case 'Faces'
        %if isempty(obj.Faces), return; end
        prop = {'Faces'}; val{1} = obj.Faces';
    case 'Gradient'
        return;% let the patch object deal with it, otherwise funny effects
%         if isempty(obj.Gradient), return; end
%         prop = {'VertexNormals'}; val{1} = obj.Gradient';
    case 'Transparency'
        alpha(obj.ph,obj.Alpha);
        return;
    case 'Material'
        prop={'AmbientStrength';'DiffuseStrength';'SpecularStrength';...
              'SpecularExponent';'SpecularColorReflectance'};
          switch obj.Material
              case 'Facial'
                  val = {0.5; 0.4; 0.1; 'default'; 'default'};
              case 'Dull'
                  val = {0.3; 0.8; 0.0; 10; 1.0};
              case 'Metal'
                  val = {0.3; 0.3; 1.0; 25; 0.5};
              case 'Shiny'
                  val = {0.3; 0.6; 0.9; 20; 1.0};
              case 'Image3D'
                  val = {0.5; 0; 0; 10000000000; 0};
              case 'Default'
                  val = {'default'; 'default'; 'default'; 'default'; 'default'};
              otherwise
                  return
          end
    case 'ColorMode'
        prop = {'FaceVertexCData'};
        switch obj.ColorMode
            case 'Texture'
               val = obj.TextureColor';
            case 'Indexed'
               val = obj.IndexedColor';
            case 'Single'
                 val = repmat(obj.SingleColor,size(get(obj.ph,'Vertices'),1),1);
%                   val = repmat(obj.SingleColor,size(obj.Vertices,2),1);
            otherwise
                 return;
        end
    case 'ViewMode'
        prop = {'EdgeColor'; 'FaceColor'; 'LineStyle';'Marker';'MarkerEdgeColor';'MarkerFaceColor'};
        switch obj.ViewMode
             case 'Solid'
                 val{1,1} = 'none'; val{2,1} = 'interp'; val{3,1} = 'none';
                 val{4,1} = 'none'; val{5,1} = 'auto';   val{6,1} = 'none';
             case 'Wireframe'
                 val{1,1} = 'interp'; val{2,1} = 'none'; val{3,1} = '-';
                 val{4,1} = 'none';   val{5,1} = 'auto'; val{6,1} = 'none';
             case 'DottedWireframe'
                 val{1,1} = 'interp'; val{2,1} = 'none'; val{3,1} = '--';
                 val{4,1} = 'none';   val{5,1} = 'auto'; val{6,1} = 'none';
             case 'Solid/Wireframe'
                 val{1,1} = [0 0 0.35]; val{2,1} = 'interp'; val{3,1} = '-';
                 val{4,1} = 'none';     val{5,1} = 'auto';   val{6,1} = 'none';
             case 'Points'
                 val{1,1} = 'none'; val{2,1} = 'none'; val{3,1} = '-';
                 val{4,1} = '.';    val{5,1} = 'flat'; val{6,1} = 'flat';
             otherwise
                 return;
        end        
    case 'LightMode'
        prop = {'FaceLighting'};
        val{1} = obj.LightMode;
    case 'Visible'
        prop = {'Visible'};
        if obj.Visible, val = {'on'}; else val = {'off'}; end
    case 'MarkerSize'
        prop = 'MarkerSize';
        val{1} = obj.MarkerSize;
    otherwise
        return;
end
pvPairs = [prop,val]';
set(obj.ph,pvPairs{:});
end
