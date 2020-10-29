function lineSearch(obj)
    switch obj.LS
        case 'none'
               obj.X = obj.X + obj.S*obj.D;% updating Parameters
               if obj.LevelUpdate, return; end% specialized eval required
               eval(obj.ObjFun);% eval current situation
        case 'backtrack'
            ArmijoBacktrack(obj);
        case 'bracketing'
            switch obj.Bracketing
                case 'fminbnd'
                    myFminbnd(obj);
                case 'wolfe'
                    myWolfe(obj);
            end
        otherwise
            error('wrong type of LineSearch');
    end
end