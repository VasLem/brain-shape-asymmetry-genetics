function outu = upsampleAnovaStats(out, reducedTemplateAdjacency, landmarksIndices)
try
toupsample = permute(...
    cat(3,out.LM.I,out.LM.IF,out.LM.IP, out.LM.permIF,...
    out.LM.D,out.LM.DF, out.LM.DP,out.LM.permDF,...
    out.LM.F,out.LM.FF,out.LM.FP,out.LM.permFF),[3,1,2]);
catch
toupsample = permute(cat(3,out.LM.I,out.LM.D,out.LM.F),[3,1,2]);
outupsampledRaw = upsampleShape3D(out.Raw.F, reducedTemplateAdjacency, landmarksIndices);
end
outupsampled= upsampleShape3D(toupsample, reducedTemplateAdjacency, landmarksIndices);
[d,r,l] = size(outupsampled);
if d == 4
    outu.LM.I = reshape(outupsampled(1,:,:),r,l);
    outu.LM.D = reshape(outupsampled(2,:,:),r,l);
    outu.LM.F = reshape(outupsampled(3,:,:),r,l);
    outu.Raw.F = outupsampledRaw;
else
outu.LM.I = reshape(outupsampled(1,:,:),r,l);
outu.LM.IF = reshape(outupsampled(2,:,:),r,l);
outu.LM.IP = reshape(outupsampled(3,:,:),r,l);
outu.LM.permIF = reshape(outupsampled(4,:,:),r,l);
outu.LM.D = reshape(outupsampled(5,:,:),r,l);
outu.LM.DF = reshape(outupsampled(6,:,:),r,l);
outu.LM.DP =reshape( outupsampled(7,:,:),r,l);
outu.LM.permDF = reshape(outupsampled(8,:,:),r,l);
outu.LM.F = reshape(outupsampled(9,:,:),r,l);
outu.LM.FF= reshape(outupsampled(10,:,:),r,l);
outu.LM.FP= reshape(outupsampled(11,:,:),r,l);
outu.LM.permFF =reshape( outupsampled(12,:,:),r,l);
outu.Total = out.Total;
end