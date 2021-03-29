function outu = upsampleAnovaStats(out, reducedTemplateAdjacency, landmarksIndices)
toupsample = permute(...
    cat(3,out.LM.I,out.LM.IF,out.LM.IP, out.LM.permIF,...
    out.LM.D,out.LM.DF, out.LM.DP,out.LM.permDF,...
    out.LM.F,out.LM.FF,out.LM.FP,out.LM.permFF),[3,1,2]);

outupsampled= upsampleShape3D(toupsample, reducedTemplateAdjacency, landmarksIndices);
outu.LM.I = squeeze(outupsampled(1,:,:));
outu.LM.IF = squeeze(outupsampled(2,:,:));
outu.LM.IP = squeeze(outupsampled(3,:,:));
outu.LM.permIF = squeeze(outupsampled(4,:,:));
outu.LM.D = squeeze(outupsampled(5,:,:));
outu.LM.DF = squeeze(outupsampled(6,:,:));
outu.LM.DP =squeeze( outupsampled(7,:,:));
outu.LM.permDF = squeeze(outupsampled(8,:,:));
outu.LM.F = squeeze(outupsampled(9,:,:));
outu.LM.FF= squeeze(outupsampled(10,:,:));
outu.LM.FP= squeeze(outupsampled(11,:,:));
outu.LM.permFF =squeeze( outupsampled(12,:,:));
outu.Total = out.Total;