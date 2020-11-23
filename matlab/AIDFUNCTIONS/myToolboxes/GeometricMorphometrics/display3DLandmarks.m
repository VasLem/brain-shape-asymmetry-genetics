function display3DLandmarks(TemplateLandmarks, Landmarks)
shape = shape3D;
if class(TemplateLandmarks) == "shape3D"
    TemplateLandmarks = TemplateLandmarks.Vertices;
end
shape.Vertices = TemplateLandmarks;
shape.VertexSize = 20;
shape.SingleColor = [0 1 0];
v = viewer(shape);
if nargin == 2
    N = size(Landmarks, 3);
    for n=1:1:N
       shape = shape3D;
       shape.Vertices = squeeze(Landmarks(:,:,n));
       shape.VertexSize = 10;
       viewer(shape,v);
    end
end
drawnow;
end