function out = reOrderRightSegmentValues(in,rend)    
         out = nan*zeros(1,length(rend.UMASK));
         out(find(rend.UMASK)) = in;
         out = out(rend.ReOrder.Ind);
         out = out(find(rend.ReOrder.UMASK));
end