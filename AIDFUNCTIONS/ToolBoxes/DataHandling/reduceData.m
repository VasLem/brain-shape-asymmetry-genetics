function redData = reduceData(Data,index)
         redData = Data;
         redData.Names = redData.Names(index);
         if isfield(redData,'Shape'), redData.Shape = redData.Shape(:,index);end
         if isfield(redData,'NormShape'), redData.NormShape = redData.NormShape(:,index);end
         if isfield(redData,'TextureColor'),redData.TextureColor = redData.TextureColor(:,index);end
         if isfield(redData,'TextureMap'),redData.TextureMap = redData.TextureMap(:,index);end
         if isfield(redData,'Properties'), redData.Properties.Values = redData.Properties.Values(:,index);end
         if isfield(redData,'Dysmorphogram'), redData.Dysmorphogram = redData.Dysmorphogram(:,index);end
         if isfield(redData,'PercDysmorph'), redData.PercDysmorph = redData.PercDysmorph(:,index);end
         if isfield(redData,'BinaryDysmorphogram'), redData.BinaryDysmorphogram = redData.BinaryDysmorphogram(:,index);end
         if isfield(redData,'BinaryPercDysmorph'), redData.BinaryPercDysmorph = redData.BinaryPercDysmorph(:,index);end
         if isfield(redData,'DistanceMap'), redData.DistanceMap = redData.DistanceMap(:,index);end
         if isfield(redData,'RMSE'), redData.RMSE = redData.RMSE(:,index);end
         if isfield(redData,'ThresholdMap'), redData.ThresholdMap = redData.ThresholdMap(:,index);end
         if isfield(redData,'PercThresh'), redData.PercThresh = redData.PercThresh(:,index);end
         if isfield(redData,'NoiseLevel'), redData.NoiseLevel = redData.NoiseLevel(:,index);end
end