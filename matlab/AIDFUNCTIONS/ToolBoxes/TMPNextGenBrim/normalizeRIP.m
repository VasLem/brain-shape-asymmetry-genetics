function [out] = normalizeRIP(rip,M,var)
     out = rip;
     for b= 1:1:var.nr2Boot % testing reference faces and rescale accordingly
        Mrip = updateRIP(var.Info{b}.MDepVar,M)';
        Prip = updateRIP(var.Info{b}.PDepVar,M)';
        tmp = rip(:,b);
        range = Prip(b)-Mrip(b);
        tmp = (tmp-Mrip(b))/range;           
        out(:,b) = tmp*var.Info{b}.Range+var.Info{b}.Mel;  
     end
end