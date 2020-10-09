function out = compareMorphs(morph1,morph2)
        % Distances (morph(1) = Scan = AA = -) (Morph(2) = Norm = BB = +)
            out.Distances = vDistances(morph1,morph2);
         % Differences
            out.Differences = vDifferences(morph1,morph2);
         % Distances along the normals
            out.NormalDistances = vNormalDistances(morph2,morph1);
         % ratio of local areas
            areas2 = localAreas(morph2);
            areas1 = localAreas(morph1);
            out.AreaRatios = -1*log(areas2./areas1);
         % curvatures and ratio of curvatures
            [Cmean2,~,Dir1,Dir2]=curvature(morph2,true);
            % getting signed curvature
            Dir3 = cross(Dir1',Dir2');
            angles = vectorAngle(morph2.Gradient,Dir3);
            signs = ones(1,length(angles));
            signs(find(angles<=90)) = -1;
            signedCmean2 = signs.*Cmean2';
            [Cmean1,~,Dir1,Dir2]=curvature(morph1,true);
            % getting signed curvature
            Dir3 = cross(Dir1',Dir2');
            angles = vectorAngle(morph1.Gradient,Dir3);
            signs = ones(1,length(angles));
            signs(find(angles<=90)) = -1;
            signedCmean1 = signs.*Cmean1';
            % smoothing of curvatures
            Cmean2 = smoothFunction(morph2,Cmean2',2,'functiondistance');
            signedCmean2 = smoothFunction(morph2,signedCmean2',2,'functiondistance');
            Cmean1 = smoothFunction(morph1,Cmean1',2,'functiondistance');
            signedCmean1 = smoothFunction(morph1,signedCmean1',2,'functiondistance');
            % computing the ratios
            out.Curvature1 = Cmean1;
            out.signedCurvature1 = signedCmean1;
            out.Curvature2 = Cmean2; 
            out.signedCurvature2 = signedCmean2;
            out.CurvatureRatios = -1*log(Cmean2./Cmean1);
            out.signedCurvatureDiff = signedCmean1-signedCmean2;% positieve increased convexity, negative is decreased convexity
end


%         % Distances (morph(1) = Scan = AA = -) (Morph(2) = Norm = BB = +)
% %             out.Distances = vDistances(morphs(1).scan,morphs(2).scan);
%          % Differences
% %             out.Differences = vDifferences(morphs(1).scan,morphs(2).scan);
%          % Distances along the normals
%             out.NormalDistances = vNormalDistances(morphs(2).scan,morphs(1).scan);
%          % ratio of local areas
% %             areas2 = localAreas(morphs(2).scan);
% %             areas1 = localAreas(morphs(1).scan);
% %             out.AreaRatios = -1*log(areas2./areas1);
%          % curvatures and ratio of curvatures
% %             Cmean2=curvature(morphs(2).scan,true);
% %             Cmean1=curvature(morphs(1).scan,true);
%             % smoothing of curvatures
% %             Cmean2 = smoothFunction(morphs(2).scan,Cmean2',2,'functiondistance');
% %             Cmean1 = smoothFunction(morphs(1).scan,Cmean1',2,'functiondistance');
%             % computing the ratios
% %             out.CurvatureRatios = -1*log(Cmean2./Cmean1);
% %             out.Curvature1 = Cmean1;
% %             out.Curvature2 = Cmean2;