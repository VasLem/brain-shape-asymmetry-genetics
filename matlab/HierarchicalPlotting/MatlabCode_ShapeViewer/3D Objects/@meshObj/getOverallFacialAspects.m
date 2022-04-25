function out = getOverallFacialAspects(obj,show)
         % aligning principal axes face with coordinate system
         obj = clone(obj);
         centerVertices(obj);
         [EigVec,~,EigVal] = my_princomp(obj.Vertices','econ');
         AxesEigVal = EigVal;
         ind = zeros(1,3);
         for i=1:1:3
             [~,ind(i)] = max(abs(EigVec(:,i)));
             AxesEigVal(ind(i)) = EigVal(i);
         end
         out.H = AxesEigVal(2);
         out.W = AxesEigVal(1);
         out.D = AxesEigVal(3);
         out.HWratio = AxesEigVal(2)/AxesEigVal(1);
         EigVec = EigVec(:,ind);
         RM = inv(EigVec);
         obj.Vertices = RM*obj.Vertices;     
         border(obj)
         BordPoints = obj.Border.Vertices(1:2,:);
         x = BordPoints(1,:);
         y = BordPoints(2,:);
         % roundess
         [xc,yc,Re,~] = circfit(x,y);
         distances = abs((x-xc).^2 + (y-yc).^2 - Re^2);
         if show
            viewer(obj); 
            figure;plot(x,y,'b.');
            lm = LMObj('Vertices',obj.Border.Vertices);
            lm.Value = distances;
            viewer(lm);
            lm.ColorMode = 'Indexed';
         end
         out.Roundness = sqrt(mean(distances.^2));
         % squarness
         bbox = [min(x) max(x) min(y) max(y)];
         distances = zeros(4,length(x));
         distances(1,:) = abs(x-bbox(1));
         distances(2,:) = abs(x-bbox(2));
         distances(3,:) = abs(y-bbox(3));
         distances(4,:) = abs(y-bbox(4));
         distances = min(distances);
         
         out.Squarness = sqrt(mean(distances.^2));
         figure;hist(distances);
         
         lm = LMObj('Vertices',obj.Border.Vertices);
         lm.Value = distances;
         viewer(lm);
         lm.ColorMode = 'Indexed';
         
end