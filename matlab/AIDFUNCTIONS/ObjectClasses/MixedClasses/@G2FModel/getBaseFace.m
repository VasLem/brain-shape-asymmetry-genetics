function getBaseFace(obj,pred)
            A = [obj.RIPS obj.RIPA];
            B = obj.Model.Tcoeff;
            Aval = [pred.RIPS pred.RIPA];
            [~,~,~,~,MC] = plsregress(A,B,size(A,2));
            mC = mean(A);
            mShape = mean(B);
            DX = mC-Aval;
            DY = DX*MC(2:end,:);
            NY = mShape-DY;
            pred.BaseFace = getScan(obj.Model,NY');
end