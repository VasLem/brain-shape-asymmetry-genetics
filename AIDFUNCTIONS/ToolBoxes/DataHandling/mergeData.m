function Data = mergeData(Data1,Data2)
         %Data = Data1;
         Data.Names = [Data1.Names Data2.Names];
         if isfield(Data1,'Shape'), Data.Shape = [Data1.Shape Data2.Shape];end
         if isfield(Data1,'NormShape'), Data.NormShape = [Data1.NormShape Data2.NormShape];end
%         if isfield(Data1,'TextureColor'),Data.TextureColor = [Data1.TextureColor Data2.TextureColor];end
%         if isfield(Data1,'TextureMap'),Data.TextureMap = [Data1.TextureMap Data2.TextureMap];end
%          if isfield(Data,'Properties'), Data.Shape = [Data.Shape Data2.Shape];end
%          if isfield(Data,'Dysmorphogram'), Data.Shape = [Data.Shape Data2.Shape];end
%          if isfield(Data,'PercDysmorph'), Data.Shape = [Data.Shape Data2.Shape];end
%          if isfield(Data,'BinaryDysmorphogram'), Data.Shape = [Data.Shape Data2.Shape];end
%          if isfield(Data,'BinaryPercDysmorph'), Data.Shape = [Data.Shape Data2.Shape];end
%          if isfield(Data,'DistanceMap'), rData.Shape = [Data.Shape Data2.Shape];end
%          if isfield(Data,'RMSE'), Data.Shape = [Data.Shape Data2.Shape];end
%          if isfield(Data,'ThresholdMap'), Data.Shape = [Data.Shape Data2.Shape];end
%          if isfield(Data,'PercThresh'), Data.Shape = [Data.Shape Data2.Shape];end
%          if isfield(Data,'NoiseLevel'), Data.Shape = [Data.Shape Data2.Shape];end
         % eliminate doubles:
         [~,I,~] = unique(Data.Names,'first');
         Data = reduceData(Data,I);
end