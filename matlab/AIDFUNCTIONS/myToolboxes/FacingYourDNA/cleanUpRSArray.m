function rsout = cleanUpRSArray(rs)
        nrs = length(rs);
        rsindex = zeros(1,nrs);
        rsfind = strfind(rs,'rs');
        index1 = find(~cellfun(@isempty,rsfind));
        rsindex(index1) = cell2mat(rsfind);        
        % If any RS make sure it is moved to the starting position
        todo = find(rsindex>1);
        if ~isempty(todo)
           redrs = rs(todo);
           redrsindex = rsindex(todo);
           redrs = cellfun(@(x,c) x(c:end),redrs,num2cell(redrsindex'),'UniformOutput',false);
           rs(todo) = redrs;
           rsindex(todo) = 1;
        end
        
        pointfind = strfind(rs,':');
        index2 = find(cellfun(@isempty,pointfind));
        rsindex(index2) = 0;
        indexk = setdiff(1:nrs,index2);
        
        pointindex = zeros(nrs,1);
        pointfind = pointfind(indexk);
        pointindex(indexk) = cellfun(@(x) x(1),pointfind)-1;
        
        todo = find(rsindex);
        if ~isempty(todo)
           redrs = rs(todo);
           redpointindex = pointindex(todo);
           redrs = cellfun(@(x,c) x(1:c),redrs,num2cell(redpointindex),'UniformOutput',false);
           rs(todo) = redrs;
        end
        rsout = rs;
end