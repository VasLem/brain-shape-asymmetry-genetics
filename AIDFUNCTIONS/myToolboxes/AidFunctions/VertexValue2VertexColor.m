function out = VertexValue2VertexColor(value,clim,cmap)

         n = length(value);
         out = nan*zeros(n,3);
         TMP = value;
         TMP(value<=clim(1)) = clim(1);
         TMP
         









end