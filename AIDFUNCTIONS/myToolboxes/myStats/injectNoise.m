function X = injectNoise(X,NL,NLsets)
         % Noise injection
         X = repmat(X,[1 1 NLsets]);
         for i=2:1:NLsets
             X(:,:,i) = X(:,:,1) + NL*randn(size(X(:,:,1)));
         end  
end