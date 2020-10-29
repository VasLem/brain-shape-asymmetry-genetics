function [morphs,STATS] = createMorphs(AR,shape,RefScan,BF,CondValues)
         [AABfor,BABfor] = eliminateNAN(AR,shape);
         [~,~,~,~,M,~,~,~] = plsregress(AABfor,BABfor,size(AABfor,2));
         Y = mean(BABfor);
         %RefScan.Vertices = reshape(Y,3,size(BABfor,2)/3);
         m = mean(AABfor);
         mv = std(AABfor);
         if nargin < 5
            CondValues = m(1:end-1); 
         end
         range(1) = BF*(m(end)-3*mv(end));
         range(2) = BF*(m(end)+3*mv(end));
         counter = 0;
         for k=1:1:2
             kval = range(k);
             counter = counter+1;
             V = [CondValues kval];
             DX = m-V;
             DY = DX*M(2:end,:);
             NY = Y-DY;
             scan = clone(RefScan);
             delete(scan.PoseLM);
             scan.Vertices = reshape(NY,3,length(NY)/3);
             morphs(counter).scan = clone(scan);
             if ~isempty(RefScan.PoseLM)
                 transferPoseLM(morphs(counter).scan,RefScan);
             end
         end
end