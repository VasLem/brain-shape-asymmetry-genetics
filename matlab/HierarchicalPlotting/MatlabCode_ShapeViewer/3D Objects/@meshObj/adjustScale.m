function adjustScale(obj)
    range = median(max(obj.Location,[],2)-min(obj.Location,[],2));
    if range>= 80, factor = 1;
    elseif range >=8, factor = 10;
    elseif range >=0.8, factor = 100;
    else
      factor = 1000;
    end
    obj.Location = obj.Location*factor;
end

