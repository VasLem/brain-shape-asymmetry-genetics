function out = signedCurvaturev2(obj)
            [Cmean1,~,Dir1,Dir2]=curvaturev2(obj,true);
            % getting signed curvature
            Dir3 = cross(Dir1',Dir2');
            % angles = vectorAngle(obj.Gradient,Dir3);
            %obj.Vertices = double(obj.Vertices);
            gradient = obj.VertexNormals;
            % angles = vectorAngle(gradient,Dir3);
            angles = vectorAngle(gradient',Dir3);
            signs = ones(1,length(angles));
            signs(find(angles<=90)) = -1;
            signedCmean1 = signs.*Cmean1';
            % smoothing of curvatures
            out = smoothFunctionv2(obj,signedCmean1',2,'functiondistance');
end