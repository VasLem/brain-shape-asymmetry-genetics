function smoothSurface(obj,naver,type,varargin)

  obj.Vertices = smoothFunction(obj,obj.Vertices,naver,type,varargin{:})';
  updateChildren(obj,'Smooth');

end