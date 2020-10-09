function [scan] = getBaseFace(A,Aval,B,Model)
            [A,B] = eliminateNAN(A,B);
            [~,~,~,~,MC] = plsregress(A,B,size(A,2));
            mC = mean(A);
            mShape = mean(B);
            DX = mC-Aval;
            DY = DX*MC(2:end,:);
            NY = mShape-DY;
            scan = getScan(Model,NY');
end