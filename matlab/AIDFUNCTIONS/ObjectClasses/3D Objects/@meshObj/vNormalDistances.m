function out = vNormalDistances(startobj,endobj)
         VF = vectorField3D;
         VF.StartPoints = startobj.Vertices;
         VF.EndPoints = endobj.Vertices;
         [~,cosangle] = vectorAngle(-1*startobj.Gradient,VF.Direction);
         out = VF.Length.*cosangle;
end