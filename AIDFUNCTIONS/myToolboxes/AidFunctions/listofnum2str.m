function list = listofnum2str(list)
    nData = length(list);
    for i=1:1:nData
       list{i} = num2str(list{i});  
    end
end