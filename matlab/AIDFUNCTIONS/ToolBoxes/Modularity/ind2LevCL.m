function [lev,cl] = ind2LevCL(ind,CLUSTERS)      
        maxind = zeros(1,length(CLUSTERS));
        for i=1:1:length(CLUSTERS)
            maxind(i) = sum(CLUSTERS(1:i));
        end
        lev = find(ind<=maxind);
        lev = lev(1);
        if lev==1; lev = 0; cl = 0; return; end
        lev = lev-1;
        cl = ind-maxind(lev);
end