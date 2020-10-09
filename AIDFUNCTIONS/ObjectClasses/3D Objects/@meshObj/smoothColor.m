function smoothColor(obj,naver,type,varargin)  
  switch obj.ColorMode
      case 'Texture'
        if strcmp(type,'distance'), type = 'texturedistance'; end  
        obj.TextureColor = smoothFunction(obj,obj.TextureColor,naver,type,varargin{:})';
      case 'Indexed'
        if strcmp(type,'distance'), type = 'indexeddistance'; end  
        obj.IndexedColor = smoothFunction(obj,obj.IndexedColor,naver,type,varargin{:})';  
      case 'Single'
          return;
  end

end