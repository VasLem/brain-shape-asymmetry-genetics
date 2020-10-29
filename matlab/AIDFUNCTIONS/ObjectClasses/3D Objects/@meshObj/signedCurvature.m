function out = signedCurvature(obj)
            [Cmean1,~,Dir1,Dir2]=curvature(obj,true);
            % getting signed curvature
            Dir3 = cross(Dir1',Dir2');
            angles = vectorAngle(obj.Gradient,Dir3);
            signs = ones(1,length(angles));
            signs(find(angles<=90)) = -1;
            signedCmean1 = signs.*Cmean1';
            % smoothing of curvatures
            out = smoothFunction(obj,signedCmean1',2,'functiondistance');
end