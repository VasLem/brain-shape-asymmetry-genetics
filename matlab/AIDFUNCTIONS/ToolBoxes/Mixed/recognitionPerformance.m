function [out,fig] = recognitionPerformance(Model,Est,B)
         n = size(Est,1);        
         rank = zeros(1,n);
         f = statusbar('RecPerf');drawnow;
         for i=1:1:n             
             dist = zeros(1,n);
             for j=1:1:n
                 [~,~,percout] = getAngle(Model,Est(i,:),B(j,:),'Mahalanobis');
                 dist(j) = 100-percout;
             end
             [~,I] = sort(dist);
             rank(i) = (find(I==i)/n)*100;
             statusbar(i/n,f);drawnow;
         end
         delete(f);
         out.rank = rank;
         X = (1:1:100);
         for i=1:1:100
             Y(i) = (sum(rank<=i)/n)*100;
         end
         out.X = X;
         out.Y = Y;
         fig = figure;plot(X,Y);
end