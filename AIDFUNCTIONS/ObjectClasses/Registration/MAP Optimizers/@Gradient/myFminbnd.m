function myFminbnd(obj)
     X = obj.X;
     D = obj.D;
     s = 0;
     opt.BracketStep = obj.S;
     opt.MultipleBracketStep = 0;
     [s1,s2] = fbracket(@ErrorFunction,s,opt,X,D,obj);
     %tmp = [s1 s2]
     [obj.S,e_next] = fminbnd(@ErrorFunction,s1,s2,opt,X,D,obj); %#ok<NASGU>
     %tmp = [e_next obj.F]
end

function error = ErrorFunction(S,X,D,obj)
        obj.X = X + S*D;
        fun(obj,'f');
        error = obj.F;
end