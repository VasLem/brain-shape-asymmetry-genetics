function out = extractNoseFeatures(feat,nose)        
% Nose angle
       LM = transferLM(nose,feat.RefNose,feat.AngleLM);
       a = LM.Vertices(:,2)-LM.Vertices(:,1);
       b = LM.Vertices(:,2)-LM.Vertices(:,3);
       out.Angle = acosd(dot(a,b)/(norm(a)*norm(b)));
% Power/Charm Line
       LM = transferLM(nose,feat.RefNose,feat.PowerCharmLM);
       Dmiddle = getLMDistance(LM.Vertices(:,4),LM.Vertices(:,3));
       Dextr1 = getLMDistance(LM.Vertices(:,2),LM.Vertices(:,2));
       Dextr2 = getLMDistance(LM.Vertices(:,5),LM.Vertices(:,6));
       Dextr = (Dextr1+Dextr2)/2;
       out.PowerCharm = -1*log(Dmiddle/Dextr);
% Overall profile
       LM = transferLM(nose,feat.RefNose,feat.ProfileLM);
       a = LM.Vertices(:,2)-LM.Vertices(:,1);
       b = LM.Vertices(:,1)-LM.Vertices(:,3);
       %out.Profile = acosd((dot(a,b))/(norm(a)*norm(b)));
       out.Profile = norm(cross(a,b))/norm(a);
% Profile dimple
       LM = transferLM(nose,feat.RefNose,feat.DimpleLM);
       a = LM.Vertices(:,2)-LM.Vertices(:,1);
       b = LM.Vertices(:,1)-LM.Vertices(:,4);
       out.Dimple = norm(cross(a,b))/norm(a);
% size tip
       curv = signedCurvature(nose);
       out.TipSize = mean(curv(feat.TipIndex));
% Constraint Points
       LM = transferLM(nose,feat.RefNose,feat.ConstrainLM);
       out.Constraints = LM.Vertices;
end

