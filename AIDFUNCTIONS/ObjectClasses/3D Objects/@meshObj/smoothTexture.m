function smoothTexture(obj,naver,type,varargin)
  
  if strcmp(type,'distance'), type = 'texturedistance'; end
  obj.TextureColor = smoothFunction(obj,obj.TextureColor,naver,type,varargin{:})';
end