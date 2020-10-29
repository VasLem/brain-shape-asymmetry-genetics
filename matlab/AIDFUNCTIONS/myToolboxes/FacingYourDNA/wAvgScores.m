function out = wAvgScores(in,w)
         out = in.*w;
         out = sum(out,1)./sum(w,1);
end