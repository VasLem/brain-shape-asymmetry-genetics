function p = getP(p)
% function to get the point list from objects or pass on given point list
       switch class(p)
              case {'meshObj' 'LMObj' 'floatingS'} 
                p = p.Vertices;
              otherwise
                if ~size(p,1)==3, p = p';end
        end
end