function out = TextureFromShape(obj,scan)
         if nargout == 1
            scan = clone(scan);
            out = scan;
         end
         shapecoeff = getCoeff(obj.Shape,scan);
         coeff = Vec2Coeff(obj,shapecoeff,(1:obj.nrSC));

end