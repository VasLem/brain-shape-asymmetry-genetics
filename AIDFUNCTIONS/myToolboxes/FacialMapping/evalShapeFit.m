function [out] = evalShapeFit(DATA,Template)
         N = length(DATA);
         Shapes = zeros(Template.nVertices,N);
         for i=1:1:N    
            origshape = DATA{i}.Shape;
            mappedshape = DATA{i}.MappedShape;
            [bar,index] = cart2baryKNN(origshape.Vertices,mappedshape.Vertices);
            projectedshape = clone(mappedshape);          
            [projectedshape.Vertices] = bary2cartKNN(origshape.Vertices,index,bar);
            
            Shapes(:,i) = sqrt(sum((mappedshape.Vertices-projectedshape.Vertices).^2,2));
         end
         out.Shapes = Shapes;
         out.Avg = mean(Shapes(:));
         out.std = std(Shapes(:));
         out.Perc1 = sum(Shapes(:)>1)/length(Shapes(:));
end