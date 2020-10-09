function Data = normalizePropertyModel(obj,prop,value,path)
         Data.Shape = zeros(obj.Average.nrV*3,obj.n);
         for i=1:1:obj.n
             disp(num2str(i));
             coeff = setPorpertyValue(obj,obj.Tcoeff(i,:),prop,value,path);
             scan = getScan(coeff);
             Data.Shape(:,i) = scan.Vertices(:);
         end
end