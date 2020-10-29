function Data = normalizeData2Property(obj,prop,value,path)
         Data.Shape = zeros(obj.Average.nrV*3,obj.n);
         Data.Properties.Names = obj.PropNames;
         Data.Properties.nr = length(obj.PropNames);
         Data.Properties.Values = zeros(Data.Properties.nr,obj.n);
         f = statusbar('Normalising');
         for i=1:1:obj.n
             %disp(num2str(i));
             coeff = setPropertyValue(obj,obj.Tcoeff(i,:),prop,value,path);
             scan = getScan(obj,coeff);
             Data.Shape(:,i) = scan.Vertices(:);
             Data.Properties.Values(:,i) = scan.UserData.Prop;
             clear coeff;
             delete(scan);
             statusbar(i/obj.n,f);
         end
         delete(f);
end