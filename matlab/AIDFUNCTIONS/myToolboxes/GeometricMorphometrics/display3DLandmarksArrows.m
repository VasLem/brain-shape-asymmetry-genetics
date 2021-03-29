function display3DLandmarksArrows(TemplateLandmarks, Landmarks)
if class(TemplateLandmarks) == "shape3D"
    TemplateLandmarks = TemplateLandmarks.Vertices;
end
if class(Landmarks) == "shape3D"
    Landmarks = Landmarks.Vertices;
end
diff = Landmarks - TemplateLandmarks;
figure();
quiver3(TemplateLandmarks(:,1),TemplateLandmarks(:,2),TemplateLandmarks(:,3), diff(:,1), diff(:,2), diff(:,3))
axis equal
end