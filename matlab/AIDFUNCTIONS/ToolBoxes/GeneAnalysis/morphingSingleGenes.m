function [morphs] = morphingSingleGenes(AR,shape,RefScan,BF)
         [AABfor,BABfor] = eliminateNAN(AR,shape);
         [~,~,~,~,M,~,~,~] = plsregress(AABfor,BABfor,3);
         Y = mean(BABfor);
         RefScan.Vertices = reshape(Y,3,size(BABfor,2)/3);
         m = mean(AABfor);
         mv = std(AABfor);
         range(1) = BF*(m(3)-3*mv(3));
         range(2) = BF*(m(3)+3*mv(3));
         counter = 0;
         for k=1:1:2
             kval = range(k);
             counter = counter+1;
             %V = [m(1:2) k*mv(3)];
             V = [m(1) -2 kval];
             %V = [m(1) m(2) kval];
             DX = m-V;
             DY = DX*M(2:end,:);
             NY = Y-DY;
             scan = clone(RefScan);
             scan.Vertices = reshape(NY,3,length(NY)/3);
             %save(['M2_' num2str(k)],'scan');
%              v = viewer(scan);
%              scan.Material = 'Dull';
%              scan.SingleColor = [0.7 0.7 0.7];
%              v.BackgroundColor = [1 1 1];
             morphs(counter).scan = clone(scan);
         end
end