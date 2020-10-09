function list = listofstr2num(list)
    nData = length(list);
    for i=1:1:nData
       list{i} = str2double(list{i});  
    end
end