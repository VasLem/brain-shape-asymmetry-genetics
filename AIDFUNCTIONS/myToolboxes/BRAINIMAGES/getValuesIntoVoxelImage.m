function im = getValuesIntoVoxelImage(val,BRAINMASK)
   im = nan*single(BRAINMASK.ImgCanvas);
   for i=1:1:BRAINMASK.nBrainVoxels
       pos = BRAINMASK.Position(:,BRAINMASK.Index(i));
       pos(1)= size(BRAINMASK.ImgCanvas,1)+1-pos(1);
       im(pos(1),pos(2),pos(3)) = val(i);
   end
end