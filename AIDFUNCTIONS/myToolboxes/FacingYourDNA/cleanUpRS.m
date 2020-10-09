function rs = cleanUpRS(rs)
        tmp = strfind(lower(rs),'rs');
        if~isempty(tmp)
           tmp2 = strfind(lower(rs),':');
           if ~isempty(tmp2)
              rs = rs(tmp(1):tmp2(1)-1);
           else
              rs = rs(tmp(1):end);
           end
        end
end