function out = regressionOfSlopes(B1,B2)
         % Note B1 and B2 are both from the same module, and thus have the
         % same dimensionality.
         nPRed = length(B1);
         W = [eye(nPRed);eye(nPRed)];
         bc = [B1;B2];
         out = regstats(bc,W,'linear','fstat');
end