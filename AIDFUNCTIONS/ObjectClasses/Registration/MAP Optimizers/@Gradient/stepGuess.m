function stepGuess(obj)
% Estimate the initial step for Line Search
    if obj.MstepEval == 0% First Mstep iteration
       obj.S = obj.BestStep;% Subclass dependent
       obj.Old.GtD = obj.GtD;
    else
        switch obj.LSinit
           case 'newton'
               obj.S = 1;
           case 'similar'
               obj.S = obj.S*min(2,(obj.Old.GtD)/(obj.GtD));
               obj.Old.GtD = obj.GtD;
           case 'quadratic'
                obj.S = min(1,2*(obj.F-obj.Old.F)/(obj.GtD));
                obj.Old.F = obj.F;
        end
        if obj.S <= 0
           obj.S = 1;
        end
    end
    obj.S = obj.BoundStep;% Subclass dependent
end