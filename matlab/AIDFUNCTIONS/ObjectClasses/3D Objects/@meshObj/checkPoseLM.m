function checkPoseLM(obj,nr)
       if nargin < 2, nr = 5; end
       if isempty(obj.PoseLM), indicatePoseLM(obj); end
       if ~(obj.PoseLM.nrV == nr), indicatePoseLM(obj); end
end