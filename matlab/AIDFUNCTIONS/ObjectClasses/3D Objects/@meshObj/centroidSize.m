function out = centroidSize(obj)
         gravity = mean(obj.Vertices,2);
         Differences = repmat(gravity,1,obj.nrV)-obj.Vertices;
         Distances = sqrt(sum(Differences.^2));
         out = sqrt(sum(Distances.^2));
end