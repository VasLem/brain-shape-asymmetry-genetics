function [out,f] = showShapeEffect(obj,coeff,type,name)
% showShapeEffect(obj,coeff,type)
step = 1;
  if isscalar(coeff)
     pc = coeff;
     coeff = obj.AvgCoeff;
     coeff(pc) = 1;
     step = 3*sqrt(obj.EigVal(pc));
  end
  if nargin<4, name = 'Shape Effect';end
  if nargin < 3,type = 'distance';end
  out = clone(obj.Average);
  vec = obj.EigVec*coeff*step;
  diff = Vec2Struc(obj,vec);
  VF = [];
  switch lower(type)
      case 'distance';
          out.Value = sqrt(sum(diff.^2));
      case 'x'
          out.Value = diff(1,:);
      case 'y'
          out.Value = diff(2,:);
      case 'z'
          out.Value = diff(3,:);
      case 'vector'
          startscan = clone(obj.Average);
          endscan = clone(obj.Average);
          endscan.Vertices = endscan.Vertices - 5*diff;
          showVectorField(startscan,endscan);
          return;
  end
  out.ColorMode = 'Indexed';
  out.Material = 'Dull';
  f = viewer(out);
  set(f.Figure,'Name',name);
  set(f.Toolbar.light_toggle,'State','on');
  set(f.Toolbar.link_toggle,'State','on');
  if ~isempty(VF),viewer(VF,'Viewer',f);end
end